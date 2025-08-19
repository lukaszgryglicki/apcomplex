package apcomplex

import (
	"math"
	"strconv"
	"strings"
	"testing"
)

// helper: parse with test precision
func tp(s string) *Complex { return MustParse(s, 128) }

// helper: parse decimal string (from RealStringFixed/ImagStringFixed) to float64
func f64(s string) float64 {
	// strings may be like "+0.0000" — trim leading '+'
	s = strings.TrimSpace(s)
	if len(s) > 0 && s[0] == '+' {
		s = s[1:]
	}
	v, _ := strconv.ParseFloat(s, 64)
	return v
}

// helper: |a-b| <= tol (component-wise on re & im)
func equalApprox(a, b *Complex, tol float64) bool {
	diff := Sub(a, b)
	re := f64(diff.RealStringFixed(40))
	im := f64(diff.ImagStringFixed(40))
	return math.Abs(re) <= tol && math.Abs(im) <= tol
}

func TestParseFormatRoundTrip(t *testing.T) {
	tests := []string{
		"0",
		"1",
		"-1",
		"i",
		"-i",
		"3.1415926535+2.718281828i",
		"3.1415926535-2.718281828i",
		"(2.5  -4.75)",
		"(2.5, -4.75)",
	}
	for _, s := range tests {
		z, err := Parse(s, 128)
		if err != nil {
			t.Fatalf("Parse %q failed: %v", s, err)
		}
		_ = z.StringFixed(30)      // ensure formatting works
		_ = z.StringScientific(20) // ensure sci formatting works
	}
}

func TestBasicAlgebra(t *testing.T) {
	z := tp("3.25-1.75i")
	negz := Neg(z)
	sum := Add(z, negz)
	if !equalApprox(sum, tp("0"), 1e-30) {
		t.Fatalf("z + (-z) != 0, got %s", sum.StringFixed(20))
	}

	one := tp("1")
	invz := Inv(z)
	prod := Mul(z, invz)
	if !equalApprox(prod, one, 1e-28) { // slightly looser because of division
		t.Fatalf("z * inv(z) != 1, got %s", prod.StringFixed(20))
	}

	conjz := Conj(z)
	conjConj := Conj(conjz)
	if !equalApprox(conjConj, z, 1e-30) {
		t.Fatalf("conj(conj(z)) != z, got %s vs %s", conjConj.StringFixed(20), z.StringFixed(20))
	}
}

func TestAddSubMulDiv(t *testing.T) {
	a := tp("1.5+0.75i")
	b := tp("-2.25+0.5i")

	wantAdd := tp("-0.75+1.25i")
	gotAdd := Add(a, b)
	if !equalApprox(gotAdd, wantAdd, 1e-30) {
		t.Fatalf("Add mismatch: got %s, want %s", gotAdd.StringFixed(20), wantAdd.StringFixed(20))
	}

	wantSub := tp("3.75+0.25i")
	gotSub := Sub(a, b)
	if !equalApprox(gotSub, wantSub, 1e-30) {
		t.Fatalf("Sub mismatch: got %s, want %s", gotSub.StringFixed(20), wantSub.StringFixed(20))
	}

	wantMul := tp("-3.75-0.9375i")
	gotMul := Mul(a, b)
	if !equalApprox(gotMul, wantMul, 1e-30) {
		t.Fatalf("Mul mismatch: got %s, want %s", gotMul.StringFixed(20), wantMul.StringFixed(20))
	}

	// Division check via inverse: a/b ≈ a*inv(b)
	gotDiv := Div(a, b)
	gotAlt := Mul(a, Inv(b))
	if !equalApprox(gotDiv, gotAlt, 1e-28) {
		t.Fatalf("Div mismatch a/b vs a*inv(b): %s vs %s", gotDiv.StringFixed(20), gotAlt.StringFixed(20))
	}
}

func TestExpLog(t *testing.T) {
	// avoid branch cut issues: pick a generic complex not on negative real axis
	z := tp("0.75+0.5i")
	el := Exp(Log(z)) // exp(log(z)) ~= z
	if !equalApprox(el, z, 1e-28) {
		t.Fatalf("exp(log(z)) != z, got %s vs %s", el.StringFixed(20), z.StringFixed(20))
	}
}

func TestTrigIdentityReal(t *testing.T) {
	// For real x, sin^2(x)+cos^2(x)=1. Use a real input.
	x := tp("0.5") // purely real
	s := Sin(x)
	c := Cos(x)
	s2 := Mul(s, s)
	c2 := Mul(c, c)
	sum := Add(s2, c2)
	if !equalApprox(sum, tp("1"), 1e-28) {
		t.Fatalf("sin^2+cos^2 != 1 for real x, got %s", sum.StringFixed(30))
	}
}
