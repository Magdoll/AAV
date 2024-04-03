# Build Docker images
all: laava laava_dev

clean:
	rm -r workflow-outputs/*
	rm .nextflow.log*
	rm test/build/*

laava laava_dev: %: %.dockerfile laava.conda_env.yml
	docker build -t ghcr.io/formbio/$@:latest -f $< .
