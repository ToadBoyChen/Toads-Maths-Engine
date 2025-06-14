# Use the latest Rust image as the base image
FROM rust:latest as builder

# Set the working directory inside the container
WORKDIR /usr/src/app

# Copy the Cargo.toml and Cargo.lock files first to leverage Docker caching
COPY Cargo.toml Cargo.lock ./

# Create an empty project to cache dependencies
RUN mkdir src && echo "fn main() {}" > src/main.rs && cargo build --release

# Copy the actual source code
COPY . .

# Build the application
RUN cargo build --release

# Use a minimal base image for the final container
FROM debian:buster-slim

# Set the working directory inside the container
WORKDIR /usr/src/app

# Copy the compiled binary from the builder stage
COPY --from=builder /usr/src/app/target/release/docker-rust-bin .

# Expose the port your application will run on (if applicable)
EXPOSE 8080

# Command to run the application
CMD ["./main"]