package main

import (
	"fmt"
	"hf/utils"
	"math"
	"testing"

	"gonum.org/v1/gonum/mat"
)

func TestInitialFockMatrix(t *testing.T) {
	t.Run("initial Fock", func(t *testing.T) {
		F, err := GenerateInitialFock()
		if err != nil {
			t.Error(err)
		}

		want := mat.NewDense(7, 7, []float64{
			-32.2545866, -2.7914909, -0.0000000, 0.0086110, 0.0000000, -0.1812967, -0.1812967,
			-2.7914909, -8.2368891, -0.0000000, -0.2282926, 0.0000000, -0.3857987, -0.3857987,
			-0.0000000, -0.0000000, -7.5428890, -0.0000000, 0.0000000, -0.1132121, 0.1132121,
			0.0086110, -0.2282926, -0.0000000, -7.4570295, -0.0000000, -0.1102196, -0.1102196,
			0.0000000, 0.0000000, 0.0000000, -0.0000000, -7.3471449, 0.0000000, 0.0000000,
			-0.1812967, -0.3857987, -0.1132121, -0.1102196, 0.0000000, -4.0329547, -0.0446466,
			-0.1812967, -0.3857987, 0.1132121, -0.1102196, 0.0000000, -0.0446466, -4.0329547})

		if !mat.EqualApprox(F, want, 0.001) {
			t.Errorf("Wrong S matrix, wanted \n%1.3f\n\n, got \n%1.3f\n\n", mat.Formatted(want), mat.Formatted(F))
		}
	})
}

func TestCreateCMatrix(t *testing.T) {
	t.Run("test conversion util", func(t *testing.T) {
		f, err := GenerateInitialFock()
		if err != nil {
			t.Error(err)
		}

		res := utils.ConvertDenseToSym(f)

		// TODO: Fix this test.
		fmt.Printf("Got \n%1.3f\n\n", mat.Formatted(res))
	})

	t.Run("C matrix", func(t *testing.T) {
		f, _ := GenerateInitialFock()

		C, err := CreateCMatrix(f)
		if err != nil {
			t.Error(err)
		}

		want := mat.NewDense(7, 7, []float64{
			-1.0015436, 0.2336245, 0.0000000, 0.0856842, 0.0000000, -0.0482226, -0.0000000,
			0.0071893, -1.0579388, -0.0000000, -0.3601105, -0.0000000, 0.4631213, 0.0000000,
			-0.0000000, -0.0000000, 1.0610702, 0.0000000, -0.0000000, -0.0000000, 0.2965071,
			0.0002671, -0.4272843, -0.0000000, 0.9399425, 0.0000000, 0.2129401, 0.0000000,
			0.0000000, 0.0000000, -0.0000000, 0.0000000, -1.0000000, -0.0000000, -0.0000000,
			-0.0018213, 0.1492533, -0.1377210, -0.0378579, 0.0000000, -0.7807003, -0.8501403,
			-0.0018213, 0.1492533, 0.1377210, -0.0378579, -0.0000000, -0.7807003, 0.8501403,
		})

		fmt.Printf("Wanted \n%1.3f\n\n Got \n%1.3f\n\n", mat.Formatted(want), mat.Formatted(C))
	})
}

func TestDensityMatrix(t *testing.T) {
	t.Run("density mat", func(t *testing.T) {
		f, _ := GenerateInitialFock()

		C, err := CreateCMatrix(f)
		if err != nil {
			t.Fatal(err)
		}

		D, err := DensityMatrix(C)
		if err != nil {
			t.Fatal(err)
		}

		want := mat.NewDense(7, 7, []float64{
			1.0650117, -0.2852166, -0.0000000, -0.0195534, -0.0000000, 0.0334496, 0.0334496,
			-0.2852166, 1.2489657, 0.0000000, 0.1135594, 0.0000000, -0.1442809, -0.1442809,
			-0.0000000, 0.0000000, 1.1258701, -0.0000000, -0.0000000, -0.1461317, 0.1461317,
			-0.0195534, 0.1135594, -0.0000000, 1.0660638, 0.0000000, -0.0993583, -0.0993583,
			-0.0000000, 0.0000000, -0.0000000, 0.0000000, 1.0000000, -0.0000000, -0.0000000,
			0.0334496, -0.1442809, -0.1461317, -0.0993583, -0.0000000, 0.0426802, 0.0047460,
			0.0334496, -0.1442809, 0.1461317, -0.0993583, -0.0000000, 0.0047460, 0.0426802,
		})

		if !mat.EqualApprox(D, want, 0.001) {
			t.Errorf("Wrong D matrix, wanted \n%1.3f\n\n, got \n%1.3f\n\n", mat.Formatted(want), mat.Formatted(D))
		}
	})
}

func TestGenerateHFEnergy(t *testing.T) {
	t.Run("hf energy", func(t *testing.T) {
		got, etot, err := GenerateHFEnergy()
		if err != nil {
			t.Fatal(err)
		}

		etot_want := -117.8397
		want := -125.842077437699

		if math.Abs(want-got) > 0.001 {
			t.Errorf("got %v, want %v\n", got, want)
		}

		if math.Abs(etot-etot_want) > 0.001 {
			t.Errorf("got %v, want %v\n", etot, etot_want)
		}
	})

}

func TestNewFockMatrix(t *testing.T) {
	t.Run("new fock matrix", func(t *testing.T) {
		H, err := GenerateCoreHamiltonian()
		if err != nil {
			t.Fatal(err)
		}

		f, err := GenerateInitialFock()

		C, err := CreateCMatrix(f)
		if err != nil {
			t.Fatal(err)
		}

		D, err := DensityMatrix(C)
		if err != nil {
			t.Fatal(err)
		}

		m, err := readFile(TwoElectronPath, false)
		if err != nil {
			t.Fatal(err)
		}

		TEI := TwoElectronIntegralTrans(m.(matrix))

		FPrime := NewFockMatrix(H, D, TEI)

		want := mat.NewDense(7, 7, []float64{
			-18.8132695, -4.8726875, -0.0000000, -0.0115290, 0.0000000, -0.8067323, -0.8067323,
			-4.8726875, -1.7909029, -0.0000000, -0.1808692, 0.0000000, -0.5790557, -0.5790557,
			-0.0000000, -0.0000000, 0.1939644, 0.0000000, 0.0000000, -0.1708886, 0.1708886,
			-0.0115290, -0.1808692, 0.0000000, 0.2391247, 0.0000000, -0.1828683, -0.1828683,
			0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.3091071, 0.0000000, 0.0000000,
			-0.8067323, -0.5790557, -0.1708886, -0.1828683, 0.0000000, -0.1450338, -0.1846675,
			-0.8067323, -0.5790557, 0.1708886, -0.1828683, 0.0000000, -0.1846675, -0.1450338})

		if !mat.EqualApprox(FPrime, want, 0.001) {
			t.Errorf("Wrong FPrime matrix, wanted \n%1.3f\n\n, got \n%1.3f\n\n", mat.Formatted(want), mat.Formatted(FPrime))
		}

	})
}
