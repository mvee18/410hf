package main

import "fmt"

func CoreHamiltonian(tke, ven matrix) [7][7]float64 {
	// The lengths SHOULD and MUST be the same. So we can index off of Tke.
	H := [7][7]float64{}

	for i, v := range tke {
		for j := range v {
			H[i][j] = tke[i][j] + ven[i][j]
			H[j][i] = tke[i][j] + ven[i][j]
		}
	}

	fmt.Println(H)

	return H
}

// We need to tranform the two electrons as they are into a 1D vector, using the
// compound indices formulae.
// We know that mu >= nu, lambda >= sigma, mu*nu >= lambda*sigma.
// The order is mu, nu, lambda, sigma.
func TwoElectronIntegralTrans(m matrix) vector {
	var ij, kl float64
	var mu, nu, lambda, sigma, value float64

	intermediate := make(matrix, 0)
	// final := make(vector, 406)

	test := make(matrix, 0)

	// This creates the first ij,kl index.
	for _, v := range m {
		new_row := make(vector, 3)

		mu, nu, lambda, sigma, value = v[0]-1, v[1]-1, v[2]-1, v[3]-1, v[4]

		ij = twoDims(mu, nu)
		kl = twoDims(lambda, sigma)

		new_row[0], new_row[1], new_row[2] = ij, kl, value

		intermediate = append(intermediate, new_row)
	}

	fmt.Println(intermediate)

	for _, v := range intermediate {
		new_row := make(vector, 2)
		mu, nu, value = v[0], v[1], v[2]

		ijkl := twoDims(mu, nu)

		new_row[0], new_row[1] = ijkl, value

		test = append(test, new_row)
	}

	// 7 x 7 x 7 x 7 using M = n(n+1)/2
	final := make(vector, 406)
	for _, v := range test {
		final[int(v[0])] = v[1]
	}

	fmt.Println(intermediate)
	fmt.Println(test)
	// We can iterate over the previous array to make the final compound
	// index.

	fmt.Println(final)
	return final
}

func twoDims(i, j float64) float64 {
	if i > j {
		return i*(i+1)/2 + j
	} else {
		return j*(j+1)/2 + i
	}
}
