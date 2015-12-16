#include <GL/glut.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>

using namespace std;

#define HEIGHT 600
#define WIDTH 900
#define PIECES 10
#define PI 3.14159

int picked_index;		//Variables to control moving a point
int picked = 0;

int complete = 0;		//Variable to indicate if bezier editing is done

struct Point
	{double x,y,z;};

vector<Point> ctrl_pts; 
Point sample[PIECES + 1];

// ****************************************MESH DATA STRUCTURE*********************************//

class Mesh
{
	private:
		// Number of vertices, faces and edges
		int n_pts, n_faces, n_edges;		
		// List of vertices, each with three co-ordinates
		vector<vector<double> > points;		
		// List containing number of vertices bounding each face
		vector<int> n_face_pts;				
		// List containing index vertices bounding each face
		vector<vector<int> > faces;			

	public:
		// Function to clear list of vertices
		void clear_pts()					
			{points.clear();}
		// Function to add vertex to vertex list
		void add_point(vector<double> vec)	
			{points.push_back(vec);}
		// Function to calculate total number of vertices
		void set_npts()						
			{n_pts = points.size();}
		// Compute vertices bounding each quad and store their indices
		void compute_faces()				
		{
			n_face_pts.clear();
			faces.clear();
			vector<int> face;
			int i;
			// Go through points from initial to final rotated position of curve
			for (i = 0; i < points.size() - PIECES - 1; i++) 
			{
				// Leave out final point of each curve, as it is part of only one quad
				if ((i + 1) % (PIECES + 1) != 0)			
				{
					n_face_pts.push_back(4);
					face.push_back(i);
					face.push_back(i + 1);
					face.push_back(i + PIECES + 2);
					face.push_back(i + PIECES + 1);
					faces.push_back(face);
					face.pop_back();
					face.pop_back();
					face.pop_back();
					face.pop_back();
				}
			}
			
			for (; i < points.size() - 1; i++)		//Go through points on final and initial rotated position of curve only, 
													//connecting those as well
			{
				n_face_pts.push_back(4);
				face.push_back(i);
				face.push_back(i + 1);
				face.push_back(i + PIECES + 2 - points.size());
				face.push_back(i + PIECES + 1 - points.size());
				faces.push_back(face);
				face.pop_back();
				face.pop_back();
				face.pop_back();
				face.pop_back();
			}
		}

		void set_nfaces()				//Compute number of faces
			{n_faces = faces.size();}

		void calc_nedges()				//Compute number of edges
			{n_edges = n_pts + n_faces - 2;}

		void write_file(string name)	//write stored mesh to an OFF/ply file
		{
			ofstream outfile;
			outfile.open(name.c_str());
			// First few lines to be written in this format for a ply file
			// ply
			// format ascii 1.0
			// element vertex (no. of vertices)
			// property float x
			// property float y
			// property float z
			// element face (no. of faces)
			// property list uchar int vertex_index
			// end_header
			outfile << "ply"<< endl<< "format ascii 1.0"<<endl<<"element vertex "<<n_pts<<endl<<
			"property float x\nproperty float y\nproperty float z\nelement face "<<n_faces<<endl<<
			"property list uchar int vertex_index"<<endl<<"end_header"<<endl;
			// For off file: 		
			// outfile << n_pts << " " << n_faces << " " << n_edges << endl;
			for (int i = 0; i < points.size(); i++)
				{outfile << points[i][0] << " " << points[i][1] << " " << points[i][2] << endl;}

			for (int i = 0; i < faces.size(); i++)
			{
				outfile << n_face_pts[i] << " ";
				for (int j = 0; j < faces[i].size(); j++)
					outfile << faces[i][j] << " ";
				outfile << endl;
			}
			outfile.close();
		}

		void draw_faces()		//Display boundaries of each quad forming the surface of revolution of a Bezier curve
		{
			vector<int> face;
			int i;
			for (i = 0; i < points.size() - PIECES - 1; i++) //Go through points from initial to final rotated position of curve
			{
				if ((i + 1) % (PIECES + 1) != 0)
				{
					glColor3d(0.0, 0.0, 1.0);		//Display horizontal lines in blue colour
					glBegin(GL_LINES);
					glVertex3d(points[i][0], points[i][1], points[i][2]);
					glVertex3d(points[i + PIECES + 1][0], points[i + PIECES + 1][1], points[i + PIECES + 1][2]);
					glVertex3d(points[i+1][0], points[i+1][1], points[i+1][2]);
					glVertex3d(points[i + PIECES + 2][0], points[i + PIECES + 2][1], points[i + PIECES + 2][2]);
					glEnd();
					glColor3d(0.0, 1.0, 0.0);		//Display vertical lines in green colour
					glBegin(GL_LINES);
					glVertex3d(points[i][0], points[i][1], points[i][2]);
					glVertex3d(points[i + 1][0], points[i + 1][1], points[i + 1][2]);
					glVertex3d(points[i + PIECES + 1][0], points[i + PIECES + 1][1], points[i + PIECES + 1][2]);
					glVertex3d(points[i + PIECES + 2][0], points[i + PIECES + 2][1], points[i + PIECES + 2][2]);
					glEnd();
				}
			}
			
			for (; i < points.size() - 1; i++)		//Go through points on final and initial rotated position of curve only
			{
				glColor3d(0.0, 0.0, 1.0);
				glBegin(GL_LINES);
				glVertex3d(points[i][0], points[i][1], points[i][2]);
				glVertex3d(points[i + PIECES + 1 - points.size()][0], points[i + PIECES + 1 - points.size()][1], points[i + PIECES + 1 - points.size()][2]);
				glVertex3d(points[i + 1][0], points[i + 1][1], points[i + 1][2]);
				glVertex3d(points[i + PIECES + 2 - points.size()][0], points[i + PIECES + 2 - points.size()][1], points[i + PIECES + 2 - points.size()][2]);
				glEnd();
				glColor3d(0.0, 1.0, 0.0);
				glBegin(GL_LINES);
				glVertex3d(points[i][0], points[i][1], points[i][2]);
				glVertex3d(points[i + 1][0], points[i + 1][1], points[i + 1][2]);
				glVertex3d(points[i + PIECES + 1 - points.size()][0], points[i + PIECES + 1 - points.size()][1], points[i + PIECES + 1 - points.size()][2]);
				glVertex3d(points[i + PIECES + 2 - points.size()][0], points[i + PIECES + 2 - points.size()][1], points[i + PIECES + 2 - points.size()][2]);
				glEnd();
			}
			
		}
};

Mesh mesh1;

// ********************************************PRIMITIVE DRAWING ROUTINES*********************************//

void drawCircle(double xc, double yc, double radius)
{
	//Draw a circle around each control point by approximating it with a polygon
	glPushMatrix();

	glTranslated(xc, yc, 0.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glBegin(GL_POLYGON);
	int sides = (2 * PI * radius) / 0.01;
	for (double i = 0; i < 2 * PI; i += PI / sides)
		glVertex3d(cos(i) * radius, sin(i) * radius, 0.0);
	glEnd();
	glPopMatrix();
}

void drawLine(Point a, Point b)
{
	//Draw a line between the given points
	glColor3f(0.0, 1.0, 0.0);
	glBegin(GL_LINES);
	glVertex3d(a.x, a.y, a.z);
	glVertex3d(b.x, b.y, b.z);
	glEnd();
}

// *******************************************BEZIER COMPUTATIONS*********************************//

Point Lerp(Point a, Point b, double t)
{
	//Perform linear interpolation between two points according to the given parameter t
	Point c;
	c.x = (1 - t)*a.x + t*b.x;
	c.y = (1 - t)*a.y + t*b.y;
	c.z = (1 - t)*a.z + t*b.z;
	return c;
}

void drawBezier()
{
	//Draw line between each pair of points stored in sample
	for (int i = 0; i < PIECES; i++)
		drawLine(sample[i], sample[i + 1]);
}

Point deCastlejau(double t, vector<Point> Points)
{
	if(Points.size()==1)
		return Points[0];
	else
	{
		vector<Point> newV;
		for(int i=0;i<Points.size()-1;i++)
			newV.push_back(Lerp(Points[i],Points[i+1],t));
		return deCastlejau(t,newV);
	}
}

void computeBezier()
{
	//Compute points on Bezier curve to form PIECES - no. of pieces of the curve 
	double step = 1.0 / PIECES;
	int i = 0;

	//For each t, perform recursive linear interpolation till you get a single point and save it
	for (double t = 0; t <= 1 && i <= PIECES; t += step, i++)
		{sample[i] = deCastlejau(t,ctrl_pts);}
}

// **************************MAIN DISPLAY ROUTINE WITH 2D/3D RENDERING CASES*********************************//

// Control variables for point number, curve formation, point deletion
int count = 0;
int form = 0;
int del = 0;

void display()
{
	glClearColor(0.0, 0.0, 0.0, 0.0);
	if (!complete)	//If editing is still being carried out
	{
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(0.0, 900.0, 0.0, 600.0, -900.0, 900.0);
		glMatrixMode(GL_MODELVIEW);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glLoadIdentity();

		//Plot all the control points
		for (int i = 0; i < ::count; i++)
		{
			if(picked_index==i && picked)
				glColor3f(0.0, 0.7, 0.0);
			else
				glColor3f(1.0, 0.0, 0.0);
			
			drawCircle(ctrl_pts[i].x, ctrl_pts[i].y, 8);
			glColor3f(0.0, 0.7, 0.0);
			glBegin(GL_POINTS);
			glVertex2d(ctrl_pts[i].x, ctrl_pts[i].y);
			glEnd();
		}

		//If form key has been pressed, form and draw the bezier curve
		if (form)
		{
			computeBezier();
			drawBezier();
		}
	}

	else			//If surface flag is set, the curve is now complete and can be revolved
	{
		//Change to 3D view
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(-900.0, 900.0, -600.0, 600.0, -450.0, 900.0);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		gluLookAt(0, 7.5, 25, 0, 0, 0, 0, 1, 0);
		

		vector<double> temp_pt;
		mesh1.clear_pts();
		for (int i = 0; i < 360; i += 5)		//Rotate the curve by 5 degrees each time
		{
			glPushMatrix();
			glRotated(i, 0.0, 1.0, 0.0);
			float mv[16];
			glGetFloatv(GL_MODELVIEW_MATRIX, mv);	//Get transformation matrix
			for (int j = 0; j < PIECES + 1; j++)
			{
				float xp = mv[0] * sample[j].x + mv[4] * sample[j].y + mv[8] * sample[j].z + mv[12];
				float yp = mv[1] * sample[j].x + mv[5] * sample[j].y + mv[9] * sample[j].z + mv[13];
				float zp = mv[2] * sample[j].x + mv[6] * sample[j].y + mv[10] * sample[j].z + mv[14];
				float wp = mv[3] * sample[j].x + mv[7] * sample[j].y + mv[11] * sample[j].z + mv[15];

				xp /= wp;
				yp /= wp;
				zp /= wp;
				temp_pt.push_back(xp);
				temp_pt.push_back(yp);
				temp_pt.push_back(zp);
				mesh1.add_point(temp_pt);	//Store transformed point in a mesh
				temp_pt.pop_back();
				temp_pt.pop_back();
				temp_pt.pop_back();
			}
			glPopMatrix();
		}
		mesh1.set_npts();				//Calculate no. of vertices
		mesh1.compute_faces();			//Find & save vertex combinations for each face
		mesh1.set_nfaces();				//Calculate no. of faces
		mesh1.calc_nedges();			//Calculate no. of edges
		mesh1.draw_faces();				//Display each face of the surface
		mesh1.write_file("Bezier_Mesh.ply");	//Save mesh to OFF/ply file
	}
	glFlush();
}

// **************************EVENT PROCESSING KEYS FOR USER CONTROL - MOUSE+KEYBOARD*********************************//
void processNormalKeys(unsigned char key, int xx, int yy)
{
	if(key=='F' || key=='f')
		form ^= 1;
	if(key==27)
		exit(0);
	if(key=='D' || key=='d')
		del = 1;
}

void mouse(int button, int state, int x, int y)
{
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)		//Left-click pressed
	{
		//If already present point is clicked, save index to update it on release
		for (int i = 0; i < ::count; i++)
		{
			if ((pow((double)(ctrl_pts[i].x - x), 2)) + (pow((double)(ctrl_pts[i].y - HEIGHT + y), 2)) <= pow((double)10, 2))
			{
				picked = 1;
				picked_index = i;
			}
		}

		// Save the point if not picked => creating a point
		if (!picked)
		{
				Point newp;
				newp.x = x; newp.y = HEIGHT-y, newp.z = 0;
				ctrl_pts.push_back(newp);
				::count++;
		}

		// Picked point selected and delete button pressed
		if(del && picked)
		{
			del = 0;
			ctrl_pts.erase(ctrl_pts.begin()+picked_index);
			::count--;
			picked = 0;
			picked_index = -1;
			glutPostRedisplay();
		}
	}

	if (button == GLUT_LEFT_BUTTON && state == GLUT_UP)			//Left-release
	{
		//If flag is set, a previously selected point has been moved. Update x, h-y in picked_index
		if (picked)
		{
			ctrl_pts[picked_index].x = x,ctrl_pts[picked_index].y = HEIGHT - y;
			picked = 0;
		}
	}
	
	if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN)		//Right-click
	{
		complete = 1;		//Bezier editing is over. Display surface of revolution
	}
	glutPostRedisplay();	//Redraw frame with new 3D settings
}

void activeMotion(int x, int y)
{
	//If flag is set, a previously selected point is being moved. Update x, h-y in picked_index
	if (picked)
		ctrl_pts[picked_index].x = x, ctrl_pts[picked_index].y = HEIGHT - y;
	glutPostRedisplay();
}

// **************************MAIN CALLER*********************************//

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB); // set the display mode
	glutInitWindowSize(WIDTH, HEIGHT); // set window size
	glutInitWindowPosition(50, 50);
	glutCreateWindow("Bezier");
	glViewport(0, 0, 900, 600);

	// register the callback functions
	glutMouseFunc(mouse);
	glutDisplayFunc(display);
	glutIdleFunc(display);
	glutMotionFunc(activeMotion);
	glutKeyboardFunc(processNormalKeys);
	glutIgnoreKeyRepeat(1);;

	// Start the interface
	glutMainLoop(); 
	return 0;
}
