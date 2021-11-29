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
	vnn, err := readFile(NuclearRepulsionEnergy, true)
	if err != nil {
		panic(err)
	}

	sint, err := readFile(OverlapIntegrals, true)
	if err != nil {
		panic(err)
	}

	ke, err := readFile(OneElectronKinetic, true)
	if err != nil {
		panic(err)
	}

	ven, err := readFile(NuclearAttraction, true)
	if err != nil {
		panic(err)
	}

	fmt.Printf("%v\n %v\n %v\n %v\n", vnn, sint, ke, ven)

	CoreHamiltonian(ke.(matrix), ven.(matrix))

	twoE, err := readFile(TwoElectronPath, false)
	if err != nil {
		panic(err)
	}

	fmt.Println(twoE)
}

// Brent says use symmetric if symmetric to put them in the right order.
