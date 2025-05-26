// This C++ file contains a very fast NMF and NNLS implementation
//
// Author:  Zach DeBruine (zacharydebruine@gmail.com)
// Source code largely derived from RcppML (github.com/zdebruine/RcppML)
// 
// Subject to terms of GPL >=2 license

#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#ifndef EIGEN_INITIALIZE_MATRICES_BY_ZERO
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
#endif

//[[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

//[[Rcpp::plugins(openmp)]]
#ifdef _OPENMP
#include <omp.h>
#endif

// SPARSE MATRIX CLASS
class spmat {
   public:
    // public member objects
    Rcpp::NumericVector x;
    Rcpp::IntegerVector i, p, Dim;

    // constructors
    spmat(Rcpp::NumericVector x, Rcpp::IntegerVector i, Rcpp::IntegerVector p, Rcpp::IntegerVector Dim) : x(x), i(i), p(p), Dim(Dim) {}
    spmat(const Rcpp::S4& s) {
        if (!s.hasSlot("x") || !s.hasSlot("i") || !s.hasSlot("p") || !s.hasSlot("Dim"))
            Rcpp::stop("provided object could not be converted to a sparse matrix in C++. Sparse matrices must generally be a Matrix::dgCMatrix");

        x = s.slot("x");
        i = s.slot("i");
        p = s.slot("p");
        Dim = s.slot("Dim");
    }
    spmat() {}

    size_t rows() { return Dim[0]; }
    size_t cols() { return Dim[1]; }

    // const column iterator
    class InnerIterator {
       public:
        InnerIterator(spmat& ptr, int col) : ptr(ptr), col_(col), index(ptr.p[col]), max_index(ptr.p[col + 1]) {}
        operator bool() const { return (index < max_index); }
        InnerIterator& operator++() {
            ++index;
            return *this;
        }
        double& value() { return ptr.x[index]; }
        int row() const { return ptr.i[index]; }

       private:
        spmat& ptr;
        int col_, index, max_index;
    };
};

// NMF HELPER FUNCTIONS
// Pearson correlation between two matrices
inline double cor(Eigen::MatrixXd& x, Eigen::MatrixXd& y) {
    double x_i, y_i, sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0, sum_y2 = 0;
    const size_t n = x.size();
    for (size_t i = 0; i < n; ++i) {
        x_i = (*(x.data() + i));
        y_i = (*(y.data() + i));
        sum_x += x_i;
        sum_y += y_i;
        sum_xy += x_i * y_i;
        sum_x2 += x_i * x_i;
        sum_y2 += y_i * y_i;
    }
    return std::abs(1 - (n * sum_xy - sum_x * sum_y) / std::sqrt((n * sum_x2 - sum_x * sum_x) * (n * sum_y2 - sum_y * sum_y)));
}

// fast symmetric matrix multiplication, A * A.transpose() - double
Eigen::MatrixXd AAt(const Eigen::MatrixXd& A) {
    Eigen::MatrixXd AAt = Eigen::MatrixXd::Zero(A.rows(), A.rows());
    AAt.selfadjointView<Eigen::Lower>().rankUpdate(A);
    AAt.triangularView<Eigen::Upper>() = AAt.transpose();
    AAt.diagonal().array() += 1e-15;
    return AAt;
}

// scale rows in w (or h) to sum to 1 and put previous rowsums in d
void scale(Eigen::MatrixXd& w, Eigen::VectorXd& d) {
    d = w.rowwise().sum();
    d.array() += 1e-15;
    for (size_t i = 0; i < w.rows(); ++i)
        for (size_t j = 0; j < w.cols(); ++j)
            w(i, j) /= d(i);
};

// calculate sort index of vector "d" in decreasing order
inline std::vector<int> sort_index(const Eigen::VectorXd& d) {
    std::vector<int> idx(d.size());
    std::iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(), [&d](size_t i1, size_t i2) { return d[i1] > d[i2]; });
    return idx;
}

// reorder rows in dynamic matrix "x" by integer vector "ind"
inline Eigen::MatrixXd reorder_rows(const Eigen::MatrixXd& x, const std::vector<int>& ind) {
    Eigen::MatrixXd x_reordered(x.rows(), x.cols());
    for (unsigned int i = 0; i < ind.size(); ++i)
        x_reordered.row(i) = x.row(ind[i]);
    return x_reordered;
}

// reorder elements in vector "x" by integer vector "ind"
inline Eigen::VectorXd reorder(const Eigen::VectorXd& x, const std::vector<int>& ind) {
    Eigen::VectorXd x_reordered(x.size());
    for (unsigned int i = 0; i < ind.size(); ++i)
        x_reordered(i) = x(ind[i]);
    return x_reordered;
}

// NNLS SOLVER
// optimized and modified from github.com/linxihui/NNLM "c_nnls" function
inline void nnls(Eigen::MatrixXd& a, Eigen::VectorXd& b, Eigen::MatrixXd& h, const size_t sample) {
    double tol = 1;
    for (uint8_t it = 0; it < 100 && (tol / b.size()) > 1e-8; ++it) {
        tol = 0;
        for (size_t i = 0; i < h.rows(); ++i) {
            double diff = b(i) / a(i, i);
            if (-diff > h(i, sample)) {
                if (h(i, sample) != 0) {
                    b -= a.col(i) * -h(i, sample);
                    tol = 1;
                    h(i, sample) = 0;
                }
            } else if (diff != 0) {
                h(i, sample) += diff;
                b -= a.col(i) * diff;
                tol += std::abs(diff / (h(i, sample) + 1e-15));
            }
        }
    }
}

// NMF PROJECTION
void c_predict(spmat A, const Eigen::MatrixXd& w, Eigen::MatrixXd& h, const double L1, const double L2, const int threads) {
    Eigen::MatrixXd a = AAt(w);
    a.diagonal().array() *= (1 - L2);
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads)
#endif
    for (size_t i = 0; i < h.cols(); ++i) {
        if (A.p[i] == A.p[i + 1]) continue;
        Eigen::VectorXd b = Eigen::VectorXd::Zero(h.rows());
        for (spmat::InnerIterator it(A, i); it; ++it)
            b += it.value() * w.col(it.row());
        b.array() -= L1;
        nnls(a, b, h, i);
    }
}

//[[Rcpp::export(rng = FALSE)]]
Eigen::MatrixXd predict_nmf(Rcpp::S4& A_, Eigen::MatrixXd& w, const double L1, const double L2, const int threads) {
    spmat A(A_);
    Eigen::MatrixXd h(w.rows(), A.cols());
    if (w.rows() == A.rows() && w.cols() != A.rows())
        w = w.transpose();

    c_predict(A, w, h, L1, L2, threads);
    return h;
}

// NMF FUNCTION
//[[Rcpp::export(rng = FALSE)]]
Rcpp::List run_nmf(const Rcpp::S4& A_, const Rcpp::S4& At_, const double tol, const uint16_t maxit, const bool verbose,
                   const double L1, const double L2, const uint16_t threads, Eigen::MatrixXd w) {
    spmat A(A_), At(At_);

    // check validity of parameters
    if (L1 >= 1 || L2 >= 1 || L1 < 0 || L2 < 0)
        Rcpp::stop("L1 and L2 must be strictly in the range (0,1]");

    if (A.rows() != At.cols() || A.cols() != At.rows())
        Rcpp::stop("A and At are not transpose-identical");

    if (w.rows() == A.rows())
        w = w.transpose();
    else if (w.cols() != A.rows())
        Rcpp::stop("dimensions of A and w are incompatible!");

    if (verbose) Rprintf("\n%4s | %8s \n---------------\n", "iter", "tol");
    Eigen::MatrixXd h(w.rows(), A.cols());
    Eigen::VectorXd d(w.rows());

    double tol_ = 1;

    // alternating least squares updates of h and then w
    for (uint16_t iter_ = 0; iter_ < maxit && tol_ > tol; ++iter_) {
        Eigen::MatrixXd w_it = w;
        c_predict(A, w, h, L1, L2, threads);  // update h
        scale(h, d);
        Rcpp::checkUserInterrupt();
        c_predict(At, h, w, L1, L2, threads);  // update w
        scale(w, d);
        // calculate tolerance of the model fit to detect convergence
        tol_ = cor(w, w_it);  // absolute correlation between "w" across consecutive iterations
        if (verbose) Rprintf("%4d | %8.2e\n", iter_ + 1, tol_);
        Rcpp::checkUserInterrupt();
    }

    // sort factors in the model by diagonal weight
    std::vector<int> indx = sort_index(d);
    w = reorder_rows(w, indx);
    d = reorder(d, indx);
    h = reorder_rows(h, indx);

    return Rcpp::List::create(
        Rcpp::Named("w") = w.transpose(),
        Rcpp::Named("d") = d,
        Rcpp::Named("h") = h);
}