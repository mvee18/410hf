package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"gonum.org/v1/gonum/mat"
)

const (
	NuclearRepulsionEnergy string = "./input/h2o/STO-3G/enuc.dat"
	OverlapIntegrals       string = "./input/h2o/STO-3G/s.dat"
	OneElectronKinetic     string = "./input/h2o/STO-3G/t.dat"
	NuclearAttraction      string = "./input/h2o/STO-3G/v.dat"
	TwoElectronPath        string = "./input/h2o/STO-3G/eri.dat"
)

type matter interface {
	iterate() (int, float64)
}

type vector []float64
type matrix [][]float64

func (v vector) iterate() (int, float64) {
	for i, val := range v {
		return i, val
	}

	return 0, 0
}

func (m matrix) iterate() (int, float64) {
	for _, val := range m {
		for i2, k := range val {
			return i2, k
		}
	}

	return 0, 0
}

func (m matrix) convert() vector {
	out := make(vector, len(m))
	for _, v := range m {
		i, j, value := v[0]-1, v[1]-1, v[2]

		ij := twoDims(i, j)

		out[int(ij)] = value
	}

	return out
}

func convertStringSlice(s []string) (vector, error) {

	vect := make(vector, len(s))

	for i, v := range s {
		fl, err := strconv.ParseFloat(v, 64)
		if err != nil {
			return nil, err
		}

		vect[i] = fl
	}

	return vect, nil

}

func readFile(fp string, splitBool bool) (matter, error) {

	ap, err := filepath.Abs(fp)
	if err != nil {
		return nil, err
	}

	f, err := os.Open(ap)
	if err != nil {
		return nil, err
	}

	defer f.Close()

	scanner := bufio.NewScanner(f)

	m := make(matrix, 0)
	for scanner.Scan() {
		trimText := strings.TrimSpace(scanner.Text())

		split := strings.Fields(trimText)
		if len(split) == 1 {
			vec := make(vector, 0)
			vnn, err := strconv.ParseFloat(trimText, 64)
			if err != nil {
				panic(err)
			}

			vec = append(vec, vnn)

			return vec, nil
		}

		vec, err := convertStringSlice(split)
		if err != nil {
			return nil, err
		}

		m = append(m, vec)
	}
	if splitBool {
		return splitSlices(m), nil
	} else {
		return m, nil
	}
}

func getValues(m matrix) vector {
	vec := make(vector, 0)
	for _, v := range m {
		for j, v2 := range v {
			if j == len(v)-1 {
				vec = append(vec, v2)
			}
		}
	}

	return vec
}

func splitSlices(m matrix) matrix {
	m_new := matrix{}

	vec := make(vector, 0)
	for _, v := range m {
		for j, v2 := range v {
			if j == len(v)-1 {
				vec = append(vec, v2)
			}
		}
	}

	c1, c2 := 0, 0
	boxes := make(vector, 0)

	for x := range vec {
		c2 += 1
		boxes = append(boxes, vec[x])

		if c2 > c1 {
			m_new = append(m_new, boxes)
			boxes = make(vector, 0)
			c1 += 1
			c2 = 0
		}
	}

	return m_new
}

func main() {
	H, err := GenerateCoreHamiltonian()

	twoE, err := readFile(TwoElectronPath, false)
	if err != nil {
		panic(err)
	}

	TEI := TwoElectronIntegralTrans(twoE.(matrix))

	S, err := generateS()
	if err != nil {
		panic(err)
	}

	Fo, err := GenerateInitialFock()
	if err != nil {
		panic(err)
	}

	C, err := CreateCMatrix(Fo)
	if err != nil {
		panic(err)
	}

	D, err := DensityMatrix(C)
	if err != nil {
		panic(err)
	}

	E := calcHFEnergy(D, H)

	Etotal, err := CalcTotalEnergy(E)
	if err != nil {
		panic(err)
	}

	final, err := RunHF(H, Fo, D, TEI, Etotal, S)
	if err != nil {
		panic(err)
	}

	fmt.Printf("The Etot is %v in MAIN\n", final)
}

func RunHF(H *mat.Dense, F *mat.Dense, D *mat.Dense, TEI []float64, Estart float64, S *mat.Dense) (float64, error) {
	converged := false

	var Eprevious float64
	Eprevious = Estart
	for !converged {
		Fnew := NewFockMatrix(H, D, TEI)

		FPrime := FockMatrix(S, Fnew)

		Cnew, err := CreateCMatrix(FPrime)
		if err != nil {
			return 0.0, err
		}

		Dnew, err := DensityMatrix(Cnew)
		if err != nil {
			return 0.0, err
		}

		E := CalcEnergyIter(Dnew, FPrime, H)

		Etotal, err := CalcTotalEnergy(E)
		if err != nil {
			panic(err)
		}

		if math.Abs(ComputeDensityDifference(Dnew, D)) < 0.0000000000001 && math.Abs(ComputeElectronicDiff(Etotal, Eprevious)) < 0.0000000000001 {
			converged = true
			// fmt.Println(ComputeDensityDifference(Dnew, D), ComputeElectronicDiff(Etotal, Eprevious))
		} else {
			fmt.Printf("%.8f %.8f %.8f\n", Etotal, ComputeDensityDifference(Dnew, D), ComputeElectronicDiff(Etotal, Eprevious))
			F = FPrime
			D = Dnew
			Eprevious = Etotal
		}

	}

	return 0.0, nil
}

func ComputeDensityDifference(d1, d2 *mat.Dense) float64 {
	rows, cols := d1.Dims()

	var RMSD float64
	for mu := 0; mu < rows; mu++ {
		for nu := 0; nu < cols; nu++ {
			RMSD += math.Pow((d1.At(mu, nu) - d2.At(mu, nu)), 2)
		}
	}

	return RMSD
}

func ComputeElectronicDiff(e1, e2 float64) float64 {
	//fmt.Printf("Previous: %v New: %v\n", e2, e1)
	return e1 - e2
}

func CalcEnergyIter(d *mat.Dense, f *mat.Dense, h *mat.Dense) float64 {
	rows, cols := d.Dims()

	var E float64
	for mu := 0; mu < rows; mu++ {
		for nu := 0; nu < cols; nu++ {
			E += d.At(mu, nu) * (h.At(mu, nu) + f.At(mu, nu))
		}
	}

	return E
}
