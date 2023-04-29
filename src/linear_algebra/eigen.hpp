#pragma once

#include <Eigen/Dense>

//FIXME namespace chem {  FIXME
namespace linear_algebra {


template <typename T>
struct EigenLib
{
    using Scalar = T;

    template <size_t Rows>
    using Vector = Eigen::Matrix<T, Rows, 1>;

    template <size_t Rows, size_t Cols>
    using Matrix = Eigen::Matrix<T, Rows, Cols>;

    template <typename Matrix>
    using LUDecomp = Eigen::FullPivLU<Eigen::Ref<Matrix>>;

    // y <- x
    template <typename Vector>
    static void copy(Vector & y, const Vector & x) {
        y = x;
    }

    // y <- alpha * y
    template <typename Vector>
    static void scale(Vector & y, const Scalar alpha) {
        y *= alpha;
    }

    // y <- alpha * (y - x)
    template <typename Vector>
    static void aymx(Vector & y, const Scalar alpha, const Vector & x) {
        y = alpha * (y - x);
    }

    // y <- alpha * x + y
    template <typename Vector>
    static void axpy(Vector & y, const Scalar alpha, const Vector & x) {
        y = alpha * x + y;
    }

    // B <- I*alpha - A
    template <typename Matrix>
    static void iama(Matrix & B, Scalar const alpha, const Matrix & A) {
        B = Matrix::Identity() * alpha - A;
    }

    // Inplace LU decomposition
    template <typename Matrix>
    static LUDecomp<Matrix> lu_decomposition(Matrix & A) {
         return LUDecomp<Matrix>(A);
    }

    // Recalculate the LU decomposition
    template <typename Decomp, typename Matrix>
    static void update_decomposition(Decomp & decomp, Matrix & A) {
        decomp.compute(A);
    }

    // Returns `true` if matrix decomposition is non-singular
    template <typename Decomp>
    static bool invertible(const Decomp & decomp) {
        return decomp.isInvertible();
    }

    // In-place matrix solution
    template <typename Decomp, typename Vector>
    static void solve(Decomp & decomp, Vector & b) {
        b = decomp.solve(b);
    }
};

}