package main

import (
	"fmt"
	"log"
	// "math"
)

func s(x, θ float64) float64 {
	return ((1.0 + x) / 2.0 * θ) + ((x - 1.0) / 2.0)
}

func F(j int, x float64) (prod float64) {
	prod = 1.0
	for m := 0; m < N; m++ {
		if j == m {
			continue
		}
		prod *= (x - lob_pts[m])
		prod /= (lob_pts[j] - lob_pts[m])
	}
	return
}

func H(v *Vector, g *Vector, A *Matrix) {
	var s Vector
	v.SetInto(&s)

	A.Mul(v)
	v.Negate()
	v.Add(g)
	for i := 0; i < N; i++ {
		v[i] = psi(v[i])
	}
	v.Negate()
	v.Add(&s)
}

func main() {
	var A Matrix

	var sum float64
	for m := 0; m < N; m++ {
		for n := 0; n < N; n++ {
			sum = 0
			for j := 0; j < N; j++ {
				sv := s(lob_pts[n], lob_pts[j])
				sum += K(lob_pts[n], sv) * F(m, sv) * lob_wts[j]
			}
			A[n][m] = ((1.0 + lob_pts[n]) / 2.0) * sum
		}
	}

	var g Vector
	for i := 0; i < N; i++ {
		g[i] = gf(lob_pts[i])
	}

	var (
		v      Vector
		vp, sa Vector
		jac    Matrix
		d      float64
		eps    float64 = 1e-15
	)

	for i := 0; i < N; i++ {
		v[i] = 0
	}

	v.SetInto(&vp)
	for {
		v.SetInto(&sa)

		H(&v, &g, &A)
		v.Negate()

		for n := 0; n < N; n++ {
			for m := 0; m < N; m++ {
				sum = 0
				for k := 0; k < N; k++ {
					sum += A[n][k] * sa[k]
				}
				jac[n][m] = psid(g[n]-sum) * A[n][m]
			}
		}
		for i := 0; i < N; i++ {
			jac[i][i] += 1
		}

		jac.Invert()
		jac.Mul(&v)
		v.Add(&sa)

		d = v.Dist(&vp)
		if d > 1e10 {
			log.Fatal("Blew up")
		}
		if d < eps {
			break
		}

		v.SetInto(&vp)
	}

	//build out solution function
	soln := func(x float64) float64 {
		sum := float64(0)
		for j := 0; j < N; j++ {
			sv := s(x, lob_pts[j])
			for k := 0; k < N; k++ {
				sum += K(x, sv) * F(k, sv) * lob_wts[j] * v[k]
			}
		}
		return gf(x) - (((1.0 + x) / 2.0) * sum)
	}

	//init our mesh
	const num_pts = 10
	mesh := [num_pts + 1]float64{}
	for i := 0; i <= num_pts; i++ {
		mesh[i] = (float64(i)/num_pts)*2.0 - 1.0
	}

	for _, pt := range mesh {
		sol := soln(pt)
		fmt.Printf("% 10.8f\t", pt)
		fmt.Printf("% 10.8f\t", sol)
		fmt.Printf("% 20.20f\n", sol-ex(pt))
	}
}

