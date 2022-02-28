//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;

//' @title Distance Covariance
//'
//' @name dcov
//'
//' @details Implements the algorithm described in \doi{10.1016/j.csda.2019.01.016} which only has O(n log(n)) complexity.
//'
//' @param x numeric vector
//' @param y numeric vector
//'
//' @usage dcov(x,y)
//'
//' @examples \dontrun{
//'
//' set.seed(1)
//' x < -rnorm(1000)
//' y < -x ^ 2
//'
//' dcov(x, y)
//' dvov(x, x)
//' dvov(x, x)
//'
//' }
//' @export
//[[Rcpp::export]]
double dcov(const arma::colvec &x, const arma::colvec &y)
{
    int n = x.n_elem;

    // Sort elements by x order
    vec ys = y.elem(sort_index(x));
    vec xs = arma::sort(x);

    vec si = cumsum(xs);

    // Element wise multiplication
    vec ax = regspace(-(n - 2), 2, n) % xs + (si(si.n_elem - 1) - 2 * si);
    mat v = join_horiz(xs, ys, xs % ys);

    arma::umat idx(n, 2, fill::zeros);
    idx.col(0) = regspace<uvec>(1, n);

    mat iv1(n, 1, fill::zeros);
    mat iv2(n, 1, fill::zeros);
    mat iv3(n, 1, fill::zeros);
    mat iv4(n, 1, fill::zeros);

    int i = 1;
    int r = 1;
    int s = 2;

    int gap, k, st1, e1, st2, e2, kf;
    double idx1, idx2;
    arma::uvec idxr(idx.n_rows);
    mat csumv(v.n_rows + 1, v.n_cols);

    idxr = idx.col(r - 1);

    while (i < n)
    {
        gap = 2 * i;
        k = 0;
        idxr = idx.col(r - 1);
        csumv = join_vert(
            zeros(1, 3),
            cumsum(v.rows(idxr - 1)));
        for (uword &j : regspace<uvec>(0, gap, n - 1))
        {
            st1 = j;
            e1 = std::min(st1 + i - 1, n - 1);
            st2 = j + i;
            e2 = std::min(st2 + i - 1, n - 1);

            while ((st1 <= e1) && (st2 <= e2))
            {
                k += 1;
                idx1 = idxr(st1) - 1;
                idx2 = idxr(st2) - 1;
                if (ys(idx1) >= ys(idx2))
                {
                    idx(k - 1, s - 1) = idx1 + 1;
                    st1 = st1 + 1;
                }
                else
                {
                    idx(k - 1, s - 1) = idx2 + 1;
                    st2 = st2 + 1;
                    iv1(idx2, 0) = iv1(idx2) + e1 - st1 + 1;
                    iv2(idx2) = iv2(idx2) + (csumv(e1 + 1, 0) - csumv(st1, 0));
                    iv3(idx2) = iv3(idx2) + (csumv(e1 + 1, 1) - csumv(st1, 1));
                    iv4(idx2) = iv4(idx2) + (csumv(e1 + 1, 2) - csumv(st1, 2));
                }
            }
            if (st1 <= e1)
            {
                kf = k + e1 - st1 + 1;
                idx(span(k, kf - 1), s - 1) = idxr.subvec(span(st1, e1));
                k = kf;
            }
            else
            {
                if (st2 <= e2)
                {
                    kf = k + e2 - st2 + 1;
                    idx(span(k, kf - 1), s - 1) = idxr.subvec(span(st2, e2));
                    k = kf;
                }
            }
        }
        i = gap;
        r = 3 - r;
        s = 3 - s;
    }

    double covterm = as_scalar(n * (x - mean(x)).t() * (y - mean(x)));
    double c1 = as_scalar(iv1.t() * v.col(2));
    double c2 = accu(iv4);
    double c3 = as_scalar(iv2.t() * ys);
    double c4 = as_scalar(iv3.t() * xs);

    // d is the Frobenus inner product of the distance matrices
    double d = 4 * ((c1 + c2) - (c3 + c4)) - 2 * covterm;

    uvec sub_y = idx.submat(regspace<uvec>(n - 1, 0),
                            regspace<uvec>(r - 1, r - 1));

    sub_y -= 1; // Due to zero indexing

    vec ysorted = ys.elem(sub_y);
    si = cumsum(ysorted);

    vec by(n, fill::zeros);
    by.elem(sub_y) = regspace(-(n - 2), 2, n) % ysorted + (si(si.n_elem - 1) - 2 * si);

    // Use doubles to avoid overflow on very large vectors
    double n1 = n * n;
    double n2 = n1 * n;
    double n3 = n2 * n;

    double term1 = d;
    double term2 = as_scalar(ax.t() * by);
    double term3 = accu(ax) * accu(by);

    // covsq equals the square of the distance covariance between x and y
    double covsq = term1 / n1 - 2 * term2 / n2 + term3 / n3;

    double cov = sqrt(covsq);

    return (cov);
}

//' @title Distance Correlation
//'
//' @name dcor
//'
//' @param x numeric vector
//' @param y numeric vector
//'
//' @usage dcor(x,y)
//'
//' @examples \dontrun{
//'
//' set.seed(1)
//' x < -rnorm(1000)
//' y < -x ^ 2
//'
//' dcor(x, y) # dcor shows dependence between x and y
//' cor(x, y) # cor does not detect any depencence due to nonlinearity
//'
//' }
//' @export
//[[Rcpp::export]]
double dcor(const arma::vec &x, const arma::vec &y)
{
    // Preperation for openmp
    dmat x_y = join_horiz(x, y);
    urowvec::fixed<4> col_idx = {0, 0, 1, 1};
    drowvec::fixed<3> dcov_vals;

    // Compute dcov(x,x), dcov(x,y) and dcov(y,y) in parallel
#pragma omp parallel for num_threads(3) schedule(static, 1)
    for (int i = 0; i < 3; i++)
    {
        dcov_vals(i) = dcov(
            x_y.col(col_idx(i)),
            x_y.col(col_idx(i + 1)));
    }

    double corr = dcov_vals(1) / sqrt(dcov_vals(0) * dcov_vals(2));
    return (corr);
}