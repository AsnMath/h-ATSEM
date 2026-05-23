ifeq ($(OS), Windows_NT)
	RM = rmdir /S /Q
	MKDIR = mkdir
else
	RM = rm -rvf
	MKDIR = mkdir -p
endif

CXX				:= g++
CXX_FLAGS 		:= -O3 -fopenmp --std=c++20
BUILD    		:= ./build
OBJ_DIR  		:= $(BUILD)/objects
APP_DIR  		:= $(BUILD)/apps
TARGET   		:= main.exe
INCLUDE  		:=
LIB  			:=
SRC      		:= $(wildcard *.cpp)

OBJECTS  		:= $(SRC:%.cpp=$(OBJ_DIR)/%.o)
DEPENDENCIES	:= $(OBJECTS:.o=.d)

all: build $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: %.cpp | $(OBJ_DIR)
	@echo "Compiling $<..."
	$(CXX) -c $< -o $@ $(CXX_FLAGS) $(INCLUDE) -MMD -MP

$(APP_DIR)/$(TARGET): $(OBJECTS) | $(APP_DIR)
	@echo "Linking $@..."
	$(CXX) -o $@ $^ $(LIB) $(CXX_FLAGS)

-include $(DEPENDENCIES)

.PHONY: run all build clean debug

$(OBJ_DIR):
	-$(MKDIR) "$@"

$(APP_DIR):
	-$(MKDIR) "$@"

build: $(OBJ_DIR) $(APP_DIR)

run: all
	./$(APP_DIR)/$(TARGET)

clean:
	-$(RM) "$(BUILD)"
