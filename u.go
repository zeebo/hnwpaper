package main

import "math"

func ex(x float64) float64 {
	return math.Exp(-1*x)
}

func psi(x float64) float64 {
	return x * x
}

func psid(x float64) float64 {
	return 2.0 * x
}

func K(x1, x2 float64) float64 {
	return x1 - x2
}

func gf(x float64) float64 {
	return math.Exp(-1*x) + (.25 * math.Exp(-2*x)) + (.5 * x) - .25
}
