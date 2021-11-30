package main

import "testing"

func TestInitialFockMatrix(t *testing.T) {
	t.Run("initial Fock", func(t *testing.T) {
		S, err := generateS()
		if err != nil {
			t.Fatal(err)
		}

		H, err := GenerateCoreHamiltonian()
		if err != nil {
			t.Fatal(err)
		}

		InitialFockMatrix(S, H)
	})
}
