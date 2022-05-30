/* Note: this Google copyright notice only applies to the original file, which has large sections copy-pasted here. my changes are under CC0 (public domain).

 * Copyright 2015 Google Inc.
 *
 * Use of this source code is governed by a BSD-style license that can be
 * found in the LICENSE file.
 */

/*
The official instructions don't work well. These alternative instructions are intended to be the shortest path to get a minimal setup running.
The Linux steps were run through successfully on October 2019.
The Windows steps are known to be broken; the broken part is Step 7. The Include and Library directories should be tweaked.

This was made by copy-pasting and fixing two sources: https://github.com/google/skia/tree/master/experimental/GLFWTest and https://gist.github.com/zester/5163313
Don't bother trying these two sources; neither of them works.


step 1: install glfw (on Linux, "sudo apt install libglfw3-dev" will get you an acceptable (and outdated) version. on Visual Studio 2017, you must build glfw from source, contrary to Internet claims that glfw's VS2015 pre-compiled version works.)
step 2: follow the Setting Up section at http://commondatastorage.googleapis.com/chrome-infra-docs/flat/depot_tools/docs/html/depot_tools_tutorial.html#_setting_up
step 3: if you're in Windows, you will need a copy of bash; cmd.exe will fail in a later step. on my system, a copy of bash came with my installation of Git for Windows.
step 4: follow https://skia.org/user/download, using the "Clone the Skia repository" section only. Use Bash, even if you're on Windows. the Windows check "where python" is useful because sometimes python ends up in stupid places for stupid reasons
step 5: go to https://skia.org/user/build and look at the instructions, but don't follow them.
Move forward to either Windows step 6 or Linux step 6.


Windows step 6, Visual Studio 2017:
here is where bash is required, because cmd.exe doesn't allow single quotes, which are necessary to give the VC path. the various skia_use_foo commands are necessary to stop VS from erroring out when the headers are missing
run these two commands, replacing "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC" with your own VC directory:
gn gen out/Static --args='is_official_build=true win_vc="C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC" skia_use_libpng=false skia_use_zlib=false skia_use_libwebp=false skia_enable_pdf=false skia_use_libjpeg_turbo=false skia_use_expat=false'
ninja -C out/Static

Windows step 7. Warning, this is outdated. The massive "FOLDER1\x" was changed to not be necessary, but I haven't tested the correct steps now:
add this file to a new VS project
append "FOLDER1\skia\include\core;FOLDER1\skia\include\gpu;FOLDER1\skia\include\config;FOLDER1\skia\include\utils;FOLDER2\glfw\include'" to the VC include directories of your project, where FOLDERX represents the directories you put them in.
	you must include all 4 skia folders because the files inside skia folders assume they see the other folders.
	if you're unfamiliar with how the include directory works, it's in Project->Properties, VC++ Directories, Include Directories.
append "FOLDER1\skia\out\Static;FOLDER2\glfw\src\Debug;" to Library Directories, again replacing FOLDERX with the true location. add "opengl32.lib;skia.lib;glfw3.lib;" to Linker->Input->Additional Dependencies
Set build mode to x64.
Build! This will produce a debug mode binary.
If in the future you want a release mode binary, you will need to re-build glfw in release mode, and change the glfw library folder to FOLDER2\glfw\src\Release;

Linux step 6, Ubuntu 19.04. October 18, 2019:
Run:
sudo apt install clang libjpeg-dev libicu-dev libwebp-dev
bin/gn gen out/Static --args='is_official_build=true cc="clang" cxx="clang++"'
ninja -C out/Static

Linux step 7:
download this file as "glfw_ship.cpp", and place it in the parent folder of the "skia" directory. (this just makes "-Iskia" in the right place)
g++ -g -std=c++1z glfw_ship.cpp -lskia -ldl -lpthread -ljpeg -lfreetype -lz -lpng -lglfw -lfontconfig -lwebp -lwebpmux -lwebpdemux -lGL -Iskia -Lskia/out/Static/
./a.out


eventually, you will want color-correct spaces, and there are 5 places below (Ctrl+F "enable correct color spaces"), where you should replace/uncomment lines to enable this.
warning: color-correct spaces don't work in VMWare, because mesa doesn't support it.
*/

#include "GLFW/glfw3.h"
#define SK_GL
#include "include/gpu/GrBackendSurface.h"
#include "include/gpu/GrDirectContext.h"
#include "include/gpu/gl/GrGLInterface.h"
#include "include/core/SkCanvas.h"
#include "include/core/SkColorSpace.h"
#include "include/core/SkSurface.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"
#include "fluidSolver_2.h"
#include <unistd.h>

//uncomment the two lines below to enable correct color spaces
//#define GL_FRAMEBUFFER_SRGB 0x8DB9
//#define GL_SRGB8_ALPHA8 0x8C43

#define SIZE 300

GrDirectContext* sContext = nullptr;
SkSurface* sSurface = nullptr;
const int screenWidth = 900;

void error_callback(int error, const char* description) {
	fputs(description, stderr);
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GL_TRUE);
}

void init_skia(int w, int h) {
	auto interface = GrGLMakeNativeInterface();
	sContext = GrDirectContext::MakeGL(interface).release();

	GrGLFramebufferInfo framebufferInfo;
	framebufferInfo.fFBOID = 0; // assume default framebuffer
	// We are always using OpenGL and we use RGBA8 internal format for both RGBA and BGRA configs in OpenGL.
	//(replace line below with this one to enable correct color spaces) framebufferInfo.fFormat = GL_SRGB8_ALPHA8;
	framebufferInfo.fFormat = GL_RGBA8;

	SkColorType colorType = kRGBA_8888_SkColorType;
	GrBackendRenderTarget backendRenderTarget(w, h,
		0, // sample count
		0, // stencil bits
		framebufferInfo);

	//(replace line below with this one to enable correct color spaces) sSurface = SkSurface::MakeFromBackendRenderTarget(sContext, backendRenderTarget, kBottomLeft_GrSurfaceOrigin, colorType, SkColorSpace::MakeSRGB(), nullptr).release();
	sSurface = SkSurface::MakeFromBackendRenderTarget(sContext, backendRenderTarget, kBottomLeft_GrSurfaceOrigin, colorType, nullptr, nullptr).release();
	if (sSurface == nullptr) abort();
}

void cleanup_skia() {
	delete sSurface;
	delete sContext;
}

void init_glfw() {
	glfwSetErrorCallback(error_callback);
	if (!glfwInit()) {
		exit(EXIT_FAILURE);
	}

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	//(uncomment to enable correct color spaces) glfwWindowHint(GLFW_SRGB_CAPABLE, GL_TRUE);
	glfwWindowHint(GLFW_STENCIL_BITS, 0);
	//glfwWindowHint(GLFW_ALPHA_BITS, 0);
	glfwWindowHint(GLFW_DEPTH_BITS, 0);
}

void draw_grid(float* grid, SkCanvas* canvas, int N) {
	float w = ((float) screenWidth) / N;
	SkPaint paint({0.058823,0.3686274, 0.6117647,1});
	for(int j = 0; j < N; j++) {
		float y = (float)j;
		for(int i = 0; i < N; i++) {
/*			if(grid[index(i,j,SIZE)] > 0.001) {
				//paint.setAlphaf(grid[index(i,j,SIZE)]);
				canvas->drawRect({i*w,j*w,(i+1)*w,(j+1)*w},paint);
			}*/
			float x = (float)x;
			float d00 = grid[IX(i,j)];
			SkPaint paint({1.0f-d00,1.0f, 1.0f-d00,1.0f});
			canvas->drawRect({i*w,j*w,(i+1)*w,(j+1)*w},paint);


		}
	}
}

void draw_square(SkCanvas* canvas,SkPaint* paint,  float x, float y, float w) {
	canvas->drawRect({x*w,y*w,(x+1)*w,(y+1)*w},*paint);
}


// Input management
double x = 0,y = 0;
double xO = 0,yO = 0;
bool rightPressed = false;
bool leftPressed = false;

static void mouse_position_callback(GLFWwindow* window, double xpos, double ypos) {
	xO = x; yO = y;
	x = xpos; y = ypos;
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
		rightPressed = true;
	}
    if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_RELEASE) {
		rightPressed = false;
	}
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
		rightPressed = true;
	}
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
		rightPressed = false;
	}
}

void process_input(FluidGrid* grid) {
	if (rightPressed || leftPressed) {
		int index_x = (int)(xO/screenWidth*SIZE); 
		int index_y = (int)(yO/screenWidth*SIZE); 

		if(index_x <= 0 || index_x >= SIZE - 1 || index_y <=0  || index_y >= SIZE - 1) {
			return;
		}

		int N = SIZE;
		// Add velocities
		if(rightPressed) {
			grid->u_prev[IX(index_x, index_y)] = 1.0f * (x - xO);
			grid->v_prev[IX(index_x, index_y)] = 1.0f * (y - yO);
		} else {
			// Add density
			grid->dens_prev[IX(index_x, index_y)] = 10.0f;
		}
		xO = x;
		yO = y;
		
    	add_source(N,grid->u,grid->u_prev, grid->dt);
    	add_source(N,grid->v,grid->v_prev, grid->dt);
    	add_source(N,grid->dens,grid->dens_prev, grid->dt);
	}
} 

int main(void) {
	// Setup GLFW
	GLFWwindow* window;
    init_glfw();
	window = glfwCreateWindow(screenWidth, screenWidth, "Simple example", NULL, NULL);
	if (!window) {
		glfwTerminate();
		exit(EXIT_FAILURE);
	}
	glfwMakeContextCurrent(window);
	//(uncomment to enable correct color spaces) glEnable(GL_FRAMEBUFFER_SRGB);

	init_skia(screenWidth, screenWidth);

	glfwSwapInterval(1);
	glfwSetKeyCallback(window, key_callback);
	glfwSetMouseButtonCallback(window, mouse_button_callback);
	glfwSetCursorPosCallback(window, mouse_position_callback);

	// Draw to the surface via its SkCanvas.
	SkCanvas* canvas = sSurface->getCanvas(); // We don't manage this pointer's lifetime.
	
	//float* grid_d = (float*) calloc(SIZE*size, sizeof(float));
	float* grid_d = (float*) calloc((SIZE+2)*(SIZE+2), sizeof(float)); // 2
	solver_init(grid_d, SIZE);

	// Draw initial grid
	SkPaint paint;
	paint.setColor(SK_ColorWHITE);
	canvas->drawPaint(paint);
	draw_grid(grid_d, canvas, SIZE);

	sContext->flush();
	printf("Sum grid: %f\n", sum(grid_d, SIZE));	
	glfwSwapBuffers(window);

	// 

	while (!glfwWindowShouldClose(window)) {
		//glfwWaitEvents();

		// Draw blank screen 	
		SkPaint paint;
		paint.setColor(SK_ColorWHITE);
		canvas->drawPaint(paint);

		// Get input from mouse
		process_input(grid);

		solver_step();
		draw_grid(grid_d, canvas, 30);
		sContext->flush();
		// printf("Sum grid: %f\n", sum(grid_d, SIZE));	
		glfwSwapBuffers(window);
		//sleep(1);
	}

	cleanup_skia();
	free(grid_d);
	glfwDestroyWindow(window);
	glfwTerminate();
	exit(EXIT_SUCCESS);
}
