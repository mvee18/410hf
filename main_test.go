package main

import (
	"fmt"
	"reflect"
	"testing"
)

func AssertMatrices(t *testing.T, path string, start, final []float64) {
	res, err := readFile(path)
	if err != nil {
		t.Error(err)
	}

	mat := res.(matrix)

	fmt.Println(len(mat))

	// The first line in the overlap integrals.
	if !reflect.DeepEqual(mat[0], start) {
		t.Errorf("got %v, wanted %v\n", mat[0], start)
	}

	// The final line in the overlap integrals.
	if !reflect.DeepEqual(mat[len(mat)-1], final) {
		t.Errorf("got %v, wanted %v\n", mat[len(mat)-1], final)
	}
}

func TestGetNuclearRepulsion(t *testing.T) {
	t.Run("Get vnn", func(t *testing.T) {
		res, err := readFile(NuclearRepulsionEnergy)
		if err != nil {
			t.Error(err)
		}

		want := 8.002367061810450

		if res.(vector)[0] != want {
			t.Errorf("wanted %v, got %v\n", res, want)
		}
	})

	t.Run("Get sint", func(t *testing.T) {
		want1 := []float64{1, 1, 1}
		want2 := []float64{7, 7, 1}
		AssertMatrices(t, OverlapIntegrals, want1, want2)
	})

	t.Run("get ke", func(t *testing.T) {
		want1 := []float64{1, 1, 29.003199945539588}
		want2 := []float64{7, 7, 0.760031883566609}

		AssertMatrices(t, OneElectronKinetic, want1, want2)
	})

	t.Run("get ven", func(t *testing.T) {
		want1 := []float64{1, 1, -61.580595358149914}
		want2 := []float64{7, 7, -5.300202953839793}

		AssertMatrices(t, NuclearAttraction, want1, want2)
	})
}
