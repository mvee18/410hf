package main

import "fmt"

func CoreHamiltonian(tke, ven matrix) [7][7]float64 {
	// The lengths SHOULD and MUST be the same. So we can index off of Tke.
	H := [7][7]float64{}

	for i, v := range H {
		for j := range v {
			fmt.Printf("\n%v %v %v\n", H[i][j], tke[i][j], ven[i][j])
			H[i][j] = tke[i][j] + ven[i][j]
		}
	}

	fmt.Println(H)

	return H
}
