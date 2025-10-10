
CXX       := clang++
CXXFLAGS  := -std=c++17 -g -Wall -Wextra
TARGET    := main

SRC       := $(shell find . -name '*.cpp')
OBJ       := $(patsubst ./%,build/%,$(SRC:.cpp=.o))

all: $(TARGET)

$(TARGET): $(OBJ)
	@echo "🔗 Linking $(TARGET)..."
	$(CXX) $(OBJ) -o $(TARGET)

build/%.o: %.cpp
	@mkdir -p $(dir $@)
	@echo "🧱 Compiling $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: all
	@echo "🚀 Running $(TARGET)..."
	./$(TARGET)

clean:
	@echo "🧹 Cleaning build files..."
	rm -rf build $(TARGET)

build:
	@mkdir -p build

.PHONY: all clean run build
