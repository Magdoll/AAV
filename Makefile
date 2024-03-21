# Build Docker images
all: laava laavadev laavanf

laava laavadev laavanf: %: %.dockerfile laava.conda_env.yml
	docker build -t $@:latest -f $< .
