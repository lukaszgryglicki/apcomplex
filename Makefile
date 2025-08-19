all: build

build: apcomplex.go
	go build -v

test: apcomplex.go apcomplex_test.go
	go test -v

clean:
	go clean -cache -v

fmt:
	go fmt *.go
