.PHONY: build format-lint-fix strict-lint

build:
	cargo build --release

format-lint-fix:
	@cargo fmt --all
	@cargo clippy --fix --all --allow-dirty

strict-lint:
	@cargo clippy --all-features -- -D warnings
