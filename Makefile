build:
	cargo build --release --bin bam2ab1

bai:
	cargo build --release --bin bam2ab1
	sudo cp target/release/bam2ab1 /usr/bin/

bai2:
	cargo build --release --bin bam2ab1
	cp target/release/bam2ab1 /usr/bin/