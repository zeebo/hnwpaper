package main

type Matrix [N][N]float64

func (m *Matrix) Mul(v *Vector) {
	var s Vector
	v.SetInto(&s)
	var sum float64
	for i := 0; i < N; i++ {
		sum = 0
		for j := 0; j < N; j++ {
			sum += m[i][j] * s[j]
		}
		v[i] = sum
	}
}

func (m *Matrix) DivideRow(r int, f float64) {
	for i := 0; i < N; i++ {
		m[r][i] /= f
	}
}

func (m *Matrix) SubRow(s, d int, mu float64) {
	for i := 0; i < N; i++ {
		m[d][i] -= m[s][i] * mu
	}
}

func (m *Matrix) Invert() {
	var aug Matrix
	for i := 0; i < N; i++ {
		aug[i][i] = 1
	}

	var v float64
	for i := 0; i < N; i++ {
		v = m[i][i]
		m.DivideRow(i, v)
		aug.DivideRow(i, v)
		for j := 0; j < N; j++ {
			if j == i {
				continue
			}
			v = m[j][i]
			m.SubRow(i, j, v)
			aug.SubRow(i, j, v)
		}
	}
	aug.SetInto(m)
}

func (m *Matrix) SetInto(o *Matrix) {
	for i := 0; i < N; i++ {
		for j := 0; j < N; j++ {
			o[i][j] = m[i][j]
		}
	}
}
