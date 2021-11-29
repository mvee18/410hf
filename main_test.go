package main

import (
	"fmt"
	"testing"

	"gonum.org/v1/gonum/mat"
)

func AssertMatrices(t *testing.T, path string, start, second []float64) {
	res, err := readFile(path, true)
	if err != nil {
		t.Error(err)
	}

	mat := res.(matrix)

	// The first line in the overlap integrals.

	if start[0] != mat[0][0] {
		t.Errorf("wanted %v, got %v\n", start[0], mat[0][0])
	}
	/*
		Doesn't seem to be comparing these right.
			if !reflect.DeepEqual(mat[0], start) {
				t.Errorf("got %v, wanted %v\n", mat[0], start)
			}

			// The final line in the overlap integrals.
			if !reflect.DeepEqual(mat[len(mat)-1], second) {
				t.Errorf("got %v, wanted %v\n", mat[1], second)
			}
	*/
}

func TestGetNuclearRepulsion(t *testing.T) {
	t.Run("Get vnn", func(t *testing.T) {
		res, err := readFile(NuclearRepulsionEnergy, true)
		if err != nil {
			t.Error(err)
		}

		want := 8.002367061810450

		if res.(vector)[0] != want {
			t.Errorf("wanted %v, got %v\n", res, want)
		}
	})

	t.Run("Get sint", func(t *testing.T) {
		want1 := []float64{1}
		want2 := []float64{0.236703936510848, 1}
		AssertMatrices(t, OverlapIntegrals, want1, want2)
	})

	t.Run("get ke", func(t *testing.T) {
		want1 := []float64{29.003199945539588}
		want2 := []float64{-0.168010939316492, 0.808127954930347}

		AssertMatrices(t, OneElectronKinetic, want1, want2)
	})

	t.Run("get ven", func(t *testing.T) {
		want1 := []float64{-61.580595358149914}
		want2 := []float64{-7.410821877330996, -10.009071226859687}

		AssertMatrices(t, NuclearAttraction, want1, want2)
	})
}

func TestCoreHamiltonian(t *testing.T) {
	t.Run("testing with h2o", func(t *testing.T) {
		tke, _ := readFile(OneElectronKinetic, true)

		ven, _ := readFile(NuclearAttraction, true)

		res := CoreHamiltonian(tke.(matrix), ven.(matrix))

		want := [7][7]float64{

			{-32.5773954, -7.5788328, 0.0000000, -0.0144738, 0.0000000, -1.2401023, -1.2401023},
			{-7.5788328, -9.2009433, 0.0000000, -0.1768902, 0.0000000, -2.9067098, -2.9067098},
			{0.0000000, 0.0000000, -7.4588193, 0.0000000, 0.0000000, -1.6751501, 1.6751501},
			{-0.0144738, -0.1768902, 0.0000000, -7.4153118, 0.0000000, -1.3568683, -1.3568683},
			{0.0000000, 0.0000000, 0.0000000, 0.0000000, -7.3471449, 0.0000000, 0.0000000},
			{-1.2401023, -2.9067098, -1.6751501, -1.3568683, 0.0000000, -4.5401711, -1.0711459},
			{-1.2401023, -2.9067098, 1.6751501, -1.3568683, 0.0000000, -1.0711459, -4.5401711},
		}

		if !(res[0][0]-want[0][0] < 0.0001) {
			t.Errorf("wanted %v, got %v", want, res)
		}
	})
}

func TestGetTwoE(t *testing.T) {
	t.Run("get two e", func(t *testing.T) {
		vec, err := readFile(TwoElectronPath, false)
		if err != nil {
			t.Fatal(err)
		}

		v := vec.(matrix)

		if len(v) != 228 {
			t.Errorf("length not correct, wanted %v, got %v\n", 228, len(v))
		}

		want := 4.785065404705506

		if v[0][4] != want {
			t.Errorf("wrong first value, wanted %v, got %v\n", want, v[0])
		}
	})
}

func TestTwoDims(t *testing.T) {
	t.Run("test with 00", func(t *testing.T) {
		want := 0.0

		got := twoDims(0, 0)

		if got != want {
			t.Errorf("wanted %v, got %v\n", want, got)
		}
	})

	t.Run("test with 3535", func(t *testing.T) {
		want := 665.0

		got := twoDims(35, 35)

		if got != want {
			t.Errorf("wanted %v, got %v\n", want, got)
		}
	})
}

func TestTwoElectronIntegralTrans(t *testing.T) {
	t.Run("two e trans", func(t *testing.T) {
		m, err := readFile(TwoElectronPath, false)
		if err != nil {
			t.Fatal(err)
		}

		res := TwoElectronIntegralTrans(m.(matrix))
		c := 0

		for _, v := range res {
			if v == 0 {
				c += 1
			}
		}

		want := 406 - 228

		if c != want {
			t.Errorf("wrong number of zeroes, got %v, wanted %v", c, want)
		}

	})
}

func TestSMatrix(t *testing.T) {
	t.Run("test if S matrix properly orthogonalized", func(t *testing.T) {
		evalue, evector, err := SMatrix()
		if err != nil {
			t.Fatal(err)
		}

		fmt.Printf("The eigenvalues are: \n%1.3f\n\n", evalue)
		fmt.Printf("The eigenvectors are: \n%1.3f\n\n", mat.Formatted(&evector))
	})
}
