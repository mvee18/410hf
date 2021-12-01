package utils

import "gonum.org/v1/gonum/mat"

func ConvertDenseToSym(d *mat.Dense) *mat.SymDense {
	i, j := d.Dims()

	m := mat.NewSymDense(i, nil)

	for x := 0; x < i; x++ {
		for y := 0; y < j; y++ {
			m.SetSym(x, y, d.At(x, y))
		}
	}

	return m
}
