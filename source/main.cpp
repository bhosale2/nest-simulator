/*
 * main.cpp
 *
 *  Created on: Jun 22, 2014
 *      Author: mgazzola
 */

#include "UsualHeaders.h"
#include "ArgumentParser.h"
#include "Test.h"
#include "Nest.h"
#include "Collision_test.h"
#include "Mindlin_friction.h"

#ifdef SNAKE_VIZ
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include "GLUT/glut.h"
#else
#include <GL/gl.h>
#endif
#endif

using namespace std;

Test * test = NULL;

#ifdef SNAKE_VIZ
struct VisualSupport
{
	static void display()
	{
	}

	static void idle(void)
	{
		glClear(GL_COLOR_BUFFER_BIT);
		test->run();
		glutSwapBuffers();
	}

	static void run(int argc, const char ** argv)
	{
		static bool bSetup = false;

		if (!bSetup)
		{
			setup(argc, argv);
			bSetup = true;
		}

		glutDisplayFunc(display);
		glutIdleFunc(idle);
		glutMainLoop();
	}

	static void setup(int argc,  const char ** argv)
	{
		glutInit(&argc, const_cast<char **>(argv));
		glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
		//glutInitWindowSize(1024,1024);
		glutInitWindowSize(700,700);
		glutCreateWindow("School");
		glutDisplayFunc(display);
		//glClearColor(1,1,1,1);
		glClearColor(0,0,0,1);

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(0, 1, 0, 1, 0, 1);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
	}
};
#endif

int main(const int argc, const char** argv)
{
	MRAG::ArgumentParser parser(argc, argv);
	const string studycase = parser("-study").asString();

	if( studycase == "NEST" )
		test = new Nest(argc, argv);
	else if( studycase == "COLL_TEST" )
		test = new Coll_test(argc, argv);
	else if( studycase == "MINDLIN_FRICTION_TEST" )
		test = new Mindlin_friction(argc, argv);
	else
	{
		printf("Study case not defined!\n");
		abort();
	}


	try
	{
#ifdef SNAKE_VIZ
		VisualSupport::run(argc, argv);
#else
		test->run();
#endif

		return 0;
	}
	catch(string &exc)
	{
		cout << exc << endl;
	}
	catch(char const* exc)
	{
		cout << exc << endl;
	}
}



