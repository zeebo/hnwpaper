package main

import "testing"

func TestLagrange(t *testing.T) {
	eps := 1e-20
	for i := 0; i < N; i++ {
		for j := 0; j < N; j++ {
			switch {
			case i == j:
				if F(i, lob_pts[j])-1.0 > eps {
					t.Errorf("%d,%d failed", i, i)
				}
			case i != j:
				if F(i, lob_pts[j]) > eps {
					t.Errorf("%d,%d failed", i, j)
				}
			}
		}
	}
}
