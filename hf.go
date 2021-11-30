package main

import (
	"errors"
	"fmt"
	"math"

	"gonum.org/v1/gonum/mat"
)

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

func symmetricOrtho(evalues []float64, evectors *mat.Dense) *mat.Dense {
	// Use the equation S^(-1/2) = PDP^T
	// Since symmetric, the transpose and inverse are the same.
	var S mat.Dense
	var sint mat.Dense

	m := mat.NewDense(len(evalues), len(evalues), nil)
	for i, v := range evalues {
		inv := func(k float64) float64 {
			return math.Pow(k, -0.5)
		}

		m.Set(i, i, inv(v))
	}

	// Since m is a diagonal matrix, the reciprocal of the square root of each
	// term is the M^(-1/2)

	LT := evectors.T()
	sint.Mul(evectors, m)

	S.Mul(&sint, LT)

	return &S
}

func SMatrix() ([]float64, mat.Dense, error) {
	m, err := generateSMatrix()
	if err != nil {
		return nil, mat.Dense{}, err
	}

	evalues, evectors, err := eigenS(m)
	if err != nil {
		return nil, mat.Dense{}, err
	}

	return evalues, evectors, nil
}

func generateSMatrix() (*mat.SymDense, error) {
	data, err := readFile(OverlapIntegrals, false)
	if err != nil {
		return nil, nil
	}

	m := mat.NewSymDense(7, nil)

	// We need to initialize the symmetric matrix.
	for _, v := range data.(matrix) {
		m.SetSym(int(v[0])-1, int(v[1])-1, v[2])
	}

	return m, nil
}

func eigenS(m *mat.SymDense) ([]float64, mat.Dense, error) {
	var eigsym mat.EigenSym

	ok := eigsym.Factorize(m, true)
	if !ok {
		return []float64{}, mat.Dense{}, errors.New("symmetric decomposition failed")
	}

	evalues := eigsym.Values(nil)

	var ev mat.Dense
	eigsym.VectorsTo(&ev)

	return evalues, ev, nil
}
