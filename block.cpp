#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector blockdis(NumericVector v, NumericMatrix D) {
    int p = v.size();
    int num_segments = p / 2;
	double det = 1.0;
    NumericVector results(num_segments + 1, 0.0);

    for (int i = 0; i < num_segments; i++) {
        int row_start = i * 2;
        int row_end = row_start + 1;
		NumericMatrix inverse(2, 2);

		NumericMatrix submatrix = D(Range(row_start, row_end), Range(row_start, row_end));
        NumericVector subvector = v[Range(row_start, row_end)];
		
		double determinant = submatrix(0, 0) * submatrix(1, 1) - submatrix(0, 1) * submatrix(1, 0);
		
		det = det*determinant;
  
	    inverse(0, 0) = submatrix(1, 1) / determinant;
	    inverse(0, 1) = -submatrix(0, 1) / determinant;
	    inverse(1, 0) = -submatrix(1, 0) / determinant;
	    inverse(1, 1) = submatrix(0, 0) / determinant;
		
        for (int j = 0; j < 2; j++) {
            results[i] += subvector[j] * (inverse(j, 0) * subvector[0] + inverse(j, 1) * subvector[1]);
        }
    }

    if (p % 2 != 0) {
        int last_row = p - 1;
        double last_result = v[last_row] * v[last_row]/D(last_row, last_row);
        results.push_back(last_result);
		det = det*D(last_row, last_row);
    }
	NumericVector final_result(2, 0.0);
	final_result[0] = sum(results);
	final_result[1] = det;
	
    return final_result;
}
