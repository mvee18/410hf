package main

import (
	"fmt"
	"hf/utils"

	"gonum.org/v1/gonum/mat"
)

func GenerateInitialFock() (*mat.Dense, error) {
	S, err := generateS()
	if err != nil {
		return nil, err
	}

	H, err := GenerateCoreHamiltonian()
	if err != nil {
		return nil, err
	}

	F := InitialFockMatrix(S, H)

	return F, nil
}

func InitialFockMatrix(s *mat.Dense, h *mat.Dense) *mat.Dense {
	var fint mat.Dense
	var F mat.Dense

	// Due to my own incompetence in using the matrix library, I need to convert
	// the 7 x 7 matrix to a mat.Dense.

	sT := s.T()

	fint.Mul(sT, h)

	F.Mul(&fint, s)

	return &F
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

func CreateCMatrix() (*mat.Dense, error) {
	f, err := GenerateInitialFock()
	if err != nil {
		return nil, err
	}

	fs := utils.ConvertDenseToSym(f)

	evalues, evectors, err := eigenS(fs)
	if err != nil {
		return nil, err
	}

	m := mat.NewDense(len(evalues), len(evalues), nil)
	for i, v := range evalues {
		m.Set(i, i, v)
	}

	fmt.Printf("evalues got \n%1.3f\n\n", mat.Formatted(m))
	fmt.Printf("evectors got \n%1.3f\n\n", mat.Formatted(&evectors))

	C, err := transformCEVectors(&evectors)
	if err != nil {
		return nil, err
	}

	fmt.Printf("C got \n%1.3f\n\n", mat.Formatted(C))

	return C, nil
}

func transformCEVectors(Co *mat.Dense) (*mat.Dense, error) {
	var C mat.Dense

	S, err := generateS()
	if err != nil {
		return nil, err
	}

	C.Mul(S, Co)

	fmt.Printf("C got \n%1.3f\n\n", mat.Formatted(&C))

	return &C, nil
}

func DensityMatrix(C *mat.Dense) (*mat.Dense, error) {
	rows, cols := C.Dims()

	electron_count := 10

	D := mat.NewDense(rows, cols, nil)

	for mu := 0; mu < rows; mu++ {
		for nu := 0; nu < cols; nu++ {
			var val float64
			for e := 0; e < electron_count/2; e++ {
				val += C.At(mu, e) * C.At(nu, e)
			}
			D.Set(mu, nu, val)
		}
	}

	fmt.Printf("D got \n%1.3f\n\n", mat.Formatted(D))

	return D, nil
}
