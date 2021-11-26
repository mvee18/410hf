package main

import (
	"bufio"
	"fmt"
	"os"
	"path/filepath"
	"strconv"
	"strings"
)

const (
	NuclearRepulsionEnergy string = "./input/h2o/STO-3G/enuc.dat"
	OverlapIntegrals       string = "./input/h2o/STO-3G/s.dat"
	OneElectronKinetic     string = "./input/h2o/STO-3G/t.dat"
	NuclearAttraction      string = "./input/h2o/STO-3G/v.dat"
)

type mat interface {
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

func readFile(fp string) (mat, error) {

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

	return m, nil
}

func main() {
	vnn, err := readFile(NuclearRepulsionEnergy)
	if err != nil {
		panic(err)
	}

	sint, err := readFile(OverlapIntegrals)
	if err != nil {
		panic(err)
	}

	ke, err := readFile(OneElectronKinetic)
	if err != nil {
		panic(err)
	}

	ven, err := readFile(NuclearAttraction)
	if err != nil {
		panic(err)
	}

	fmt.Printf("%v\n %v\n %v\n %v\n", vnn, sint, ke, ven)

	CoreHamiltonian(ke.(matrix), ven.(matrix))
}
