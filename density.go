package main

import (
	"fmt"

	"gonum.org/v1/gonum/mat"
)

func InitialFockMatrix(s *mat.Dense, h *mat.Dense) *mat.Dense {
	//var fint mat.Dense
	//var F mat.Dense

	// Due to my own incompetence in using the matrix library, I need to convert
	// the 7 x 7 matrix to a mat.Dense.

	fmt.Printf("new H matrix \n%1.3f\n\n", mat.Formatted(h))

	return h
}

// You can probably make this faster by initializing the size ahead of time but eh.
func convertHMatrix(m [7][7]float64) *mat.Dense {
	// ie., make([]float64, len(m)*len(m))
	out := make([]float64, 0)

	for _, v := range m {
		for _, k := range v {
			out = append(out, k)
		}
	}

	return mat.NewDense(len(m), len(m), out)
}
