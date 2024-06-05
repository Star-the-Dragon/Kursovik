#include "Render.h"

#include <windows.h>
#include <GL\GL.h>
#include <GL\GLU.h>

#include "MyOGL.h"

#include "Camera.h"
#include "Light.h"
#include "Primitives.h"
#include <time.h>
#include <string>

bool textureMode = true;
bool lightMode = true;

#define b_n 4			//(�) � ���� ���������� �� �����������������, ��� ��� �� ��� ��� ������, �� ���������� ������
#define b_m 4			//(�) � ���� ���� �����

Vector3 points[b_n][b_m];

Vector3 *dragPoint = 0;				//(�) ��� ���������� ������ ����� ��������� �������

//����� ��� ��������� ������
class CustomCamera : public Camera
{
public:
	//��������� ������
	double camDist;
	//���� �������� ������
	double fi1, fi2;

	
	//������� ������ �� ���������
	CustomCamera()
	{
		camDist = 15;
		fi1 = 1;
		fi2 = 1;
	}

	
	//������� ������� ������, ������ �� ����� ��������, ���������� �������
	void SetUpCamera()
	{
		//�������� �� ������� ������ ������
		lookPoint.setCoords(0, 0, 0);

		pos.setCoords(camDist*cos(fi2)*cos(fi1),
			camDist*cos(fi2)*sin(fi1),
			camDist*sin(fi2));

		if (cos(fi2) <= 0)
			normal.setCoords(0, 0, -1);
		else
			normal.setCoords(0, 0, 1);

		LookAt();
	}

	void CustomCamera::LookAt()
	{
		//������� ��������� ������
		gluLookAt(pos.X(), pos.Y(), pos.Z(), lookPoint.X(), lookPoint.Y(), lookPoint.Z(), normal.X(), normal.Y(), normal.Z());
	}



}  camera;   //������� ������ ������


//����� ��� ��������� �����
class CustomLight : public Light
{
public:
	CustomLight()
	{
		//��������� ������� �����
		pos = Vector3(1, 1, 3);
	}

	
	//������ ����� � ����� ��� ���������� �����, ���������� �������
	void  DrawLightGhismo()
	{
		glDisable(GL_LIGHTING);

		
		glColor3d(0.9, 0.8, 0);
		Sphere s;
		s.pos = pos;
		s.scale = s.scale*0.08;
		s.Show();
		
		if (OpenGL::isKeyPressed('G'))
		{
			glColor3d(0, 0, 0);
			//����� �� ��������� ����� �� ����������
			glBegin(GL_LINES);
			glVertex3d(pos.X(), pos.Y(), pos.Z());
			glVertex3d(pos.X(), pos.Y(), 0);
			glEnd();

			//������ ���������
			Circle c;
			c.pos.setCoords(pos.X(), pos.Y(), 0);
			c.scale = c.scale*1.5;
			c.Show();
		}

	}

	void SetUpLight()
	{
		GLfloat amb[] = { 0.2, 0.2, 0.2, 0 };
		GLfloat dif[] = { 1.0, 1.0, 1.0, 0 };
		GLfloat spec[] = { .7, .7, .7, 0 };
		GLfloat position[] = { pos.X(), pos.Y(), pos.Z(), 1. };

		// ��������� ��������� �����
		glLightfv(GL_LIGHT0, GL_POSITION, position);
		// �������������� ����������� �����
		// ������� ��������� (���������� ����)
		glLightfv(GL_LIGHT0, GL_AMBIENT, amb);
		// ��������� ������������ �����
		glLightfv(GL_LIGHT0, GL_DIFFUSE, dif);
		// ��������� ���������� ������������ �����
		glLightfv(GL_LIGHT0, GL_SPECULAR, spec);

		glEnable(GL_LIGHT0);
	}


} light;  //������� �������� �����




//������ ���������� ����
int mouseX = 0, mouseY = 0;

void mouseEvent(OpenGL *ogl, int mX, int mY)
{
	int dx = mouseX - mX;
	int dy = mouseY - mY;
	mouseX = mX;
	mouseY = mY;

	//������ ���� ������ ��� ������� ����� ������ ����
	if (OpenGL::isKeyPressed(VK_RBUTTON))
	{
		camera.fi1 += 0.01*dx;
		camera.fi2 += -0.01*dy;
	}

	
	//������� ���� �� ���������, � ����� ��� ����
	if (OpenGL::isKeyPressed('G') && !OpenGL::isKeyPressed(VK_LBUTTON))
	{
		LPPOINT POINT = new tagPOINT();
		GetCursorPos(POINT);
		ScreenToClient(ogl->getHwnd(), POINT);
		POINT->y = ogl->getHeight() - POINT->y;

		Ray r = camera.getLookRay(POINT->x, POINT->y);

		double z = light.pos.Z();

		double k = 0, x = 0, y = 0;
		if (r.direction.Z() == 0)
			k = 0;
		else
			k = (z - r.origin.Z()) / r.direction.Z();

		x = k*r.direction.X() + r.origin.X();
		y = k*r.direction.Y() + r.origin.Y();

		light.pos = Vector3(x, y, z);
	}

	if (OpenGL::isKeyPressed('G') && OpenGL::isKeyPressed(VK_LBUTTON))
	{
		light.pos = light.pos + Vector3(0, 0, 0.02*dy);
	}

	if (dragPoint)				//(�) ������ ��, ��� ��� ���))
	{
		(*dragPoint) = (*dragPoint) + Vector3(0, 0, 0.02 * dy);
	}
	
}

void mouseWheelEvent(OpenGL *ogl, int delta)
{

	if (delta < 0 && camera.camDist <= 1)
		return;
	if (delta > 0 && camera.camDist >= 100)
		return;

	camera.camDist += 0.01*delta;

}

void keyDownEvent(OpenGL *ogl, int key)
{
	if (key == 'L')
	{
		lightMode = !lightMode;
	}

	if (key == 'T')
	{
		textureMode = !textureMode;
	}

	if (key == 'R')
	{
		camera.fi1 = 1;
		camera.fi2 = 1;
		camera.camDist = 15;

		light.pos = Vector3(1, 1, 3);
	}

	if (key == 'F')
	{
		light.pos = camera.pos;
	}

	if (key == VK_LBUTTON)			//(�) ��� ����� ��� ���������, ���� �� ����� ������������ ����� � �������!
	{
		double mvMatrix[16], prMatrix[16];
		int viewPort[4];
		double x, y, z;
		glGetDoublev(GL_MODELVIEW_MATRIX, mvMatrix);
		glGetDoublev(GL_PROJECTION_MATRIX, prMatrix);
		glGetIntegerv(GL_VIEWPORT, viewPort);

		LPPOINT POINT = new tagPOINT();
		GetCursorPos(POINT);
		ScreenToClient(ogl->getHwnd(), POINT);
		POINT->y = ogl->getHeight() - POINT->y;

		for (int i=0;i<b_n;++i)
			for (int j = 0;j < b_m;++j)
			{
				gluProject(points[i][j].X(), points[i][j].Y(), points[i][j].Z(), mvMatrix, prMatrix, viewPort, &x, &y, &z);
				if ((x - POINT->x) * (x - POINT->x) + (y - POINT->y) * (y - POINT->y) < 49)
				{
					dragPoint = &points[i][j];
					break;
				}
			}
		delete POINT;
	}

	if (key == VK_RBUTTON)
	{
		dragPoint = 0;
	}
}

void keyUpEvent(OpenGL *ogl, int key)		//(�) � ����������� ����������. ��� ��� ������� ��������� �� ���������� ������, � ���������� ���
{
	
}
			
GLuint texId;
//����������� ����� ������ ��������
void initRender(OpenGL *ogl)
{
	//��������� �������

	//4 ����� �� �������� �������
	glPixelStorei(GL_UNPACK_ALIGNMENT, 4);

	//��������� ������ ��������� �������
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	//�������� ��������
	glEnable(GL_TEXTURE_2D);
	

	//������ ����������� ���������  (R G B)																						
		RGBTRIPLE* texarray;

		//������ ��������, (������*������*4      4, ���������   ����, �� ������� ������������ �� 4 ����� �� ������� �������� - R G B A)
		char* texCharArray;
		int texW, texH;
		OpenGL::LoadBMP("texture1.bmp", &texW, &texH, &texarray);					
		OpenGL::RGBtoChar(texarray, texW, texH, &texCharArray);

		//���������� �� ��� ��������
		glGenTextures(1, &texId);
		//������ ��������, ��� ��� ����� ����������� � ���������, ����� ����������� �� ����� ��
		glBindTexture(GL_TEXTURE_2D, texId);

		//��������� �������� � �����������, � ���������� ��� ������  ��� �� �����
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texW, texH, 0, GL_RGBA, GL_UNSIGNED_BYTE, texCharArray);


		//�������� ������
		free(texCharArray);
		free(texarray);

		//������� ����
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
																		

	//������ � ���� ����������� � "������"
	ogl->mainCamera = &camera;
	ogl->mainLight = &light;

	// ������������ �������� : �� ����� ����� ����� 1
	glEnable(GL_NORMALIZE);

	// ���������� ������������� ��� �����
	glEnable(GL_LINE_SMOOTH); 


	//   ������ ��������� ���������
	//  �������� GL_LIGHT_MODEL_TWO_SIDE - 
	//                0 -  ������� � ���������� �������� ���������(�� ���������), 
	//                1 - ������� � ���������� �������������� ������� ��������       
	//                �������������� ������� � ���������� ��������� ����������.    
	//  �������� GL_LIGHT_MODEL_AMBIENT - ������ ������� ���������, 
	//                �� ��������� �� ���������
	// �� ��������� (0.2, 0.2, 0.2, 1.0)

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 0);  

	srand(time(0));   // (�) ��� �������!!!

	for(int i = 0; i<b_n;++i)
		for (int j = 0; j < b_m;++j)
		{
			points[i][j] = Vector3(10.0 * i / (b_n - 1) - 5, 10.0 * j / (b_m - 1) - 5, rand() % 5 - 2.5);
		}
}



inline double fact(int n)		//(�) ������� ����������
{
	int r = 1;
	for (int i = 2; i <= n;++i)
		r *= i;
	return r;
}

inline double Bern(int i, int n, double u)		//(�) ��������� ����������, � ����� ���� ������ �� ����
{
	return (fact(n)  / (fact(i) * fact(n - i))) * pow(u, i) * pow(1 - u, n - i);
}

inline Vector3 Bezier(Vector3* points, int n, int m, double u, double v)		//(�) ��� ������ ��������� �������
{
	Vector3 res(0, 0, 0);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0;j<m;++j)
		{
			double Bi = Bern(i, n - 1, u);
			double Bj = Bern(j, m - 1, v);
			res = res + points[i * m + j] * Bi * Bj;
		}
	}
	return res;
}

Vector3 BCurve(Vector3* points, int n, double t)		//(�) ������ �������
{
	Vector3 res(0, 0, 0);
	for (int i = 0; i < n;i++)
		res = res + points[i] * Bern(i, n - 1, t);
	return res;
}

void drawCube()			//(�) ������ �����, ������� ���������
{
	Vector3 a(0, 0, 0), b(1, 0, 0), c(1, 1, 0), d(0, 1, 0);

	auto drawPlane = [&]()			//(�) ������ ���������, ����� ��� ����������� ������ ������
	{
		glPushMatrix();
		glTranslated(-.5, -.5, 0);
		glBegin(GL_QUADS);
		glNormal3d(0, 0, 1);
		glVertex3dv(a.toArray());
		glVertex3dv(b.toArray());
		glVertex3dv(c.toArray());
		glVertex3dv(d.toArray());
		glEnd();
		glPopMatrix();
	};

	glPushMatrix();					//(�) ��� ����� ������ ��������
	glTranslated(0, 0, -0.5);
	glRotated(180, 1, 0, 0);
	drawPlane();
	glPopMatrix();

	glPushMatrix();
	glTranslated(0, 0, 0.5);
	drawPlane();
	glPopMatrix();

	glPushMatrix();
	glTranslated(0, -0.5, 0);
	glRotated(90, 1, 0, 0);
	drawPlane();
	glPopMatrix();

	glPushMatrix();
	glTranslated(0, 0.5, 0);
	glRotated(-90, 1, 0, 0);
	drawPlane();
	glPopMatrix();

	glPushMatrix();
	glTranslated(0.5, 0, 0);
	glRotated(90, 0, 1, 0);
	drawPlane();
	glPopMatrix();

	glPushMatrix();
	glTranslated(-0.5, 0, 0);
	glRotated(-90, 0, 1, 0);
	drawPlane();
	glPopMatrix();

	glDisable(GL_LIGHTING);

	glBegin(GL_LINES);										//(�) ����� �������� ��� ������
	glColor3d(1, 0, 0);
	glVertex3d(0, 0, 0);
	glVertex3d(2, 0, 0);
	glColor3d(0, 1, 0);
	glVertex3d(0, 0, 0);
	glVertex3d(0, 2, 0);
	glColor3d(0, 0, 1);
	glVertex3d(0, 0, 0);
	glVertex3d(0, 0, 2);
	glEnd();

	glEnable(GL_LIGHTING);
}

float anim_h = 0.0008;	
double anim_t = anim_h; 

void Render(OpenGL *ogl)
{
	glPointSize(10.0);    //(�) ������ �����

	glDisable(GL_TEXTURE_2D);
	glDisable(GL_LIGHTING);

	glEnable(GL_DEPTH_TEST);
	if (textureMode)
		glEnable(GL_TEXTURE_2D);

	if (lightMode)
		glEnable(GL_LIGHTING);

	
	//��������������
	//glEnable(GL_BLEND);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_DST_ALPHA);
	//glDisable(GL_BLEND);


	//��������� ���������
	GLfloat amb[] = { 0.2, 0.2, 0.1, 1. };
	GLfloat dif[] = { 0.4, 0.65, 0.5, 1. };
	GLfloat spec[] = { 0.9, 0.8, 0.3, 1. };
	GLfloat sh = 0.1f * 256;


	//�������
	glMaterialfv(GL_FRONT, GL_AMBIENT, amb);
	//��������
	glMaterialfv(GL_FRONT, GL_DIFFUSE, dif);
	//����������
	glMaterialfv(GL_FRONT, GL_SPECULAR, spec);
	//������ �����
	glMaterialf(GL_FRONT, GL_SHININESS, sh);

	//���� ���� �������, ��� ����������� (����������� ���������)
	glShadeModel(GL_SMOOTH);
	//===================================
					

	double h = 0.1;

	Vector3 old1, old2;
	Vector3 a, b, c, d;

	glBegin(GL_TRIANGLES);						//(�) ��� ��� ����������� ����������� ������� � ��������� � ����������
	for (double u = h; u <= 1; u += h)
	{
		old1 = Bezier((Vector3*)points, b_n, b_m, u, 0);
		old2 = Bezier((Vector3*)points, b_n, b_m, u-h, 0);
		for (double v = h; v <= 1;v += h)
		{
			a = Bezier((Vector3*)points, b_n, b_m, u, v);
			b = Bezier((Vector3*)points, b_n, b_m, u-h, v);
			c = old2;
			d = old1;

			glNormal3dv((b - c).vectProisvedenie(b - a).toArray());
			glTexCoord2d(u, v);
			glVertex3dv(a.toArray());
			glTexCoord2d(u-h, v);
			glVertex3dv(b.toArray());
			glTexCoord2d(u-h, v-h);
			glVertex3dv(c.toArray());
			
			glNormal3dv((d - a).vectProisvedenie(d - c).toArray());
			glTexCoord2d(u, v);
			glVertex3dv(a.toArray());
			glTexCoord2d(u, v - h);
			glVertex3dv(d.toArray());
			glTexCoord2d(u - h, v - h);
			glVertex3dv(c.toArray());
			
			/*glNormal3dv((d - c).vectProisvedenie(b - c).toArray());
			glTexCoord2d(u, v);
			glVertex3dv(a.toArray());
			glTexCoord2d(u - h, v);
			glVertex3dv(b.toArray());
			glTexCoord2d(u - h, v - h);
			glVertex3dv(c.toArray());
			glTexCoord2d(u, v - h);
			glVertex3dv(d.toArray());*/
			old1 = a;
			old2 = b;
		}
	}
	glEnd();

	glDisable(GL_LIGHTING);
	glColor3d(1, 0.7, 0);
	for (int i = 1;i < b_n;++i)					//(�) ��� ����������� ����� �����������
		for (int j = 1;j < b_m;++j)
		{
			glBegin(GL_LINE_LOOP);
			glVertex3dv(points[i][j].toArray());
			glVertex3dv(points[i-1][j].toArray());
			glVertex3dv(points[i-1][j-1].toArray());
			glVertex3dv(points[i][j-1].toArray());
			glEnd();
		}

	glDisable(GL_DEPTH_TEST);
	glBegin(GL_POINTS);						//(�) ��� ����������� ����� �����������
	glColor3d(0, 1, 1);
	for (int i = 0; i < b_n;++i)
		for (int j = 0; j < b_m; ++j)
		{
			if (dragPoint && dragPoint == &points[i][j])
				continue;
			glVertex3dv(points[i][j].toArray());
		}
	if (dragPoint)
	{
		glColor3d(1, 0, 0);
		glVertex3dv(dragPoint->toArray());
	}
	glEnd();

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);


	////������� ���	
		
	Vector3 p[4];							
	p[0] = Vector3(0, 0, 0);
	p[1] = Vector3(-15, 20, 40);
	p[2] = Vector3(25, 50, -10);
	p[3] = Vector3(-25, -6, 10);
	glBegin(GL_LINE_STRIP);					//(�) ��� ������ ������, � �������� ������� ��������� ������
	glVertex3dv(p[0].toArray());
	glVertex3dv(p[1].toArray());
	glVertex3dv(p[2].toArray());
	glVertex3dv(p[3].toArray());
	glEnd();

	glPushMatrix();
	Vector3 pos = BCurve(p, 4, anim_t);						
	Vector3 pre_pos = BCurve(p, 4, anim_t - anim_h);
	Vector3 dir = (pos - pre_pos).normolize();

	//dir = Vector3(1, 0.1, 0.1);

	Vector3 orig(1, 0, 0);									//(�) ������� ������
	Vector3 rotX(dir.X(), dir.Y(), 0);
	rotX = rotX.normolize();
	double cosU = (orig.X() * rotX.X())+ (orig.Y() * rotX.Y())+ (orig.Z() * rotX.Z());			
	Vector3 vecpr = orig.vectProisvedenie(rotX);

	double sinSign = vecpr.Z() / abs(vecpr.Z());
	double U = acos(cosU) * 180 / 3.141592 * sinSign;

	double cosZU = (0* dir.X()) + (0 * dir.Y()) + (1 * dir.Z());	
	double ZU = acos(dir.Z()) * 180.0 / M_PI - 90;

	glTranslated(pos.X(), pos.Y(), pos.Z());
	glRotated(U, 0, 0, 1);
	glRotated(ZU, 0, 1, 0);						

	drawCube();									//(�) �������� �����
	glPopMatrix();

	glDisable(GL_LIGHTING);
	glColor3d(1, 0.3, 0);

	glBegin(GL_LINE);								
	glVertex3dv(pos.toArray());
	glVertex3dv((pos+dir.normolize()*3).toArray());
	

	glEnd();

	anim_t += anim_h;									//(�) ��� ��������
	if (anim_t > 1) anim_h = -anim_h;
	if (anim_t < 0) anim_h = -anim_h;

	glBegin(GL_LINE_STRIP);								//(�) ��� �������� ������
	glDisable(GL_LIGHTING);
	glColor3d(1, 0, 1);
	for (double t = 0; t <= 1;t += 0.01)				
	{
		glVertex3dv(BCurve(p, 4, t).toArray());
	}
	glEnd();
  
	//����� ��������� ������ ������ ���� ������� - ���������������, ��� �������� =)
	char _c[250];
	
	//����� ���� �������

}