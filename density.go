package main

import (
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

	F := FockMatrix(S, H)

	return F, nil
}

func FockMatrix(s *mat.Dense, h *mat.Dense) *mat.Dense {
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

func CreateCMatrix(f *mat.Dense) (*mat.Dense, error) {
	fs := utils.ConvertDenseToSym(f)

	evalues, evectors, err := eigenS(fs)
	if err != nil {
		return nil, err
	}

	m := mat.NewDense(len(evalues), len(evalues), nil)
	for i, v := range evalues {
		m.Set(i, i, v)
	}

	C, err := transformCEVectors(&evectors)
	if err != nil {
		return nil, err
	}

	return C, nil
}

func transformCEVectors(Co *mat.Dense) (*mat.Dense, error) {
	var C mat.Dense

	S, err := generateS()
	if err != nil {
		return nil, err
	}

	C.Mul(S, Co)

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

	return D, nil
}

func GenerateHFEnergy() (float64, float64, error) {
	f, _ := GenerateInitialFock()

	C, err := CreateCMatrix(f)
	if err != nil {
		return 0.0, 0.0, err
	}

	D, err := DensityMatrix(C)
	if err != nil {
		return 0.0, 0.0, err
	}

	H, err := GenerateCoreHamiltonian()
	if err != nil {
		return 0.0, 0.0, err
	}

	E := calcHFEnergy(D, H)

	Etot, nil := CalcTotalEnergy(E)
	if err != nil {
		return 0.0, 0.0, err
	}

	return E, Etot, nil
}

func calcHFEnergy(d *mat.Dense, h *mat.Dense) float64 {
	rows, cols := d.Dims()

	var E float64
	for mu := 0; mu < rows; mu++ {
		for nu := 0; nu < cols; nu++ {
			E += 2 * d.At(mu, nu) * h.At(mu, nu)
		}
	}

	return E
}

func CalcTotalEnergy(Ee float64) (float64, error) {
	Emat, err := readFile(NuclearRepulsionEnergy, true)
	if err != nil {
		return 0.0, err
	}

	En := Emat.(vector)[0]

	return (En + Ee), nil
}

func NewFockMatrix(H *mat.Dense, D *mat.Dense, TEI []float64) *mat.Dense {
	nao, _ := H.Dims()

	F := mat.NewDense(nao, nao, nil)

	F.Copy(H)

	for i := 0; i < nao; i++ {
		for j := 0; j < nao; j++ {
			for k := 0; k < nao; k++ {
				for l := 0; l < nao; l++ {
					ij := index(i, j)
					kl := index(k, l)
					ijkl := index(ij, kl)

					ik := index(i, k)
					jl := index(j, l)
					ikjl := index(ik, jl)

					sum := D.At(k, l) * (2*TEI[ijkl] - TEI[ikjl])
					F.Set(i, j, F.At(i, j)+sum)
				}
			}
		}
	}

	return F

}

func index(i, j int) int {
	if i > j {
		return i*(i+1)/2 + j
	} else {
		return j*(j+1)/2 + i
	}
}
