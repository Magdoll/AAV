# Build Docker images
all: laava laava_dev laava_nf

laava laava_dev laava_nf: %: %.dockerfile laava.conda_env.yml
	docker build -t ghcr.io/formbio/$@:latest -f $< .
