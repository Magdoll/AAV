# Build Docker images
all: laava laava_dev

clean:
	rm -fv .nextflow.log*
	rm -fv test/build/*
	rm -rf workflow-outputs/*

laava laava_dev: %: %.dockerfile laava.conda_env.yml
	docker build -t ghcr.io/formbio/$@:latest -f $< .
