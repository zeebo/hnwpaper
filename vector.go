package main

type Vector [N]float64

func (v *Vector) Add(x *Vector) {
	for i := 0; i < N; i++ {
		v[i] += x[i]
	}
}

func (v *Vector) Sub(x *Vector) {
	for i := 0; i < N; i++ {
		v[i] -= x[i]
	}
}

func (v *Vector) Negate() {
	for i := 0; i < N; i++ {
		v[i] *= -1.0
	}
}

func (v *Vector) Dist(x *Vector) (d float64) {
	for i := 0; i < N; i++ {
		d += (v[i] - x[i]) * (v[i] - x[i])
	}
	return
}

func (v *Vector) SetInto(x *Vector) {
	for i := 0; i < N; i++ {
		x[i] = v[i]
	}
}

func (v *Vector) Eval(pt float64) (sum float64) {
	for i := 0; i < N; i++ {
		sum += F(i, pt) * v[i]
	}
	return
}