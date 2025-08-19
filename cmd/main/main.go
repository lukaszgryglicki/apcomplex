package main

import (
	"flag"
	"fmt"
	"math"
	"os"

	ap "github.com/lukaszgryglicki/apcomplex"
)

func main() {
	baseStr := flag.String("base", "2+1e-100i", "complex base, e.g. \"2+1e-100i\" or \"(2 1e-100)\"")
	expStr := flag.String("exp", "2-1e-100i", "complex height, e.g. \"2-1e-100i\" or \"(2 -1e-100)\"")
	// exp := flag.Int("exp", 2048, "integer exponent")
	prec := flag.Uint("prec", 8192, "precision in bits (both real & imag)")
	digits := flag.Int("digits", -1, "digits for output; -1 = auto from precision")
	out := flag.String("out", "sci", "output mode: sci|fixed")
	flag.Parse()

	// Parse inputs at requested precision.
	z, err := ap.Parse(*baseStr, *prec)
	if err != nil {
		fmt.Fprintln(os.Stderr, "parse base:", err)
		os.Exit(1)
	}
	// n := ap.MustParse(fmt.Sprintf("%d", *exp), *prec)
	n, err := ap.Parse(*expStr, *prec)
	if err != nil {
		fmt.Fprintln(os.Stderr, "parse height:", err)
		os.Exit(1)
	}

	// Compute z^n
	res := ap.New(*prec).Pow(z, n)

	// Choose a safe number of digits for printing from the precision.
	d := *digits
	if d < 0 {
		d = int(float64(*prec)*math.Log10(2)) - 5 // ~significant digits; leave a tiny safety margin
		if d < 1 {
			d = 1
		}
		if d > 1<<20 {
			d = 1<<20 // sanity cap
		}
	}

	fmt.Printf("z = %s\n", z.StringScientific(d))
	fmt.Printf("n = %s\n", n.StringScientific(d))
	//fmt.Printf("n = %d\n", *exp)
	fmt.Printf("precision = %d bits, print digits ≈ %d\n", *prec, d)

	switch *out {
	case "fixed":
		fmt.Printf("z^n (fixed, %d fractional digits): %s\n", d, res.StringFixed(d))
	default:
		fmt.Printf("z^n (scientific, %d significant digits): %s\n", d, res.StringScientific(d))
	}

	// Also show integer-like prints of components (0 fractional digits).
	fmt.Printf("Re(z^n) (fixed, 0 frac): %s\n", res.RealStringFixed(d))
	fmt.Printf("Im(z^n) (fixed, 0 frac): %s\n", res.ImagStringFixed(d))

	// Exact big-int demo: 2^1024 — this prints all digits exactly.
	// two := ap.MustParse("2", *prec)
	// k := ap.MustParse("1024", *prec)
	// pow2 := ap.New(*prec).Pow(two, k)
	// fmt.Printf("2^1024 (exact integer): %s\n", pow2.RealStringFixed(0))
}

