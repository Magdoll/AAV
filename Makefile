# Build Docker images
all: laava laava_dev

laava laava_dev: %: %.dockerfile laava.conda_env.yml
	docker build -t ghcr.io/formbio/$@:latest -f $< .
