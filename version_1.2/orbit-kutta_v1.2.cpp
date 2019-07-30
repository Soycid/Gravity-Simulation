/*
TODO:
Add moons and other planets for more accuracy
Add option for relativistic model?
Figure out way to graph 
 >add zoom feature
 >display all metrics in units
 >mouse functionality
 >add menus for RK4 and other features like white background and trace
 >fix year
*/

#include <iostream>
#include <cmath>
#include <SFML/Graphics.hpp>
#include <string>
using namespace std;

//config data
string runge = "rk4";
bool trace = false;


//one pixel = 300k km
//one time unit = 1 s
const int WinX = 2*1920;
const int WinY = 2*1080;

const double SunX = WinX/2;
const double SunY = WinY/2;
//Gravitaional Constant in pix^3 kg^-1 s^-2
const double G = 2.497413*pow(10,-36);
//Mass of Sun in kg
const double Ms = 1.989*pow(10,30);

//Initial Position of Sattelite in pix (x, 0)
const double Xi = 147099761/300000;
//Initial Velocity of Sattelite in pix/s
const double Vi = pow(G*Ms*(2/Xi - 1/((1.496/3)*pow(10,3))),.5);
//Step size in seconds
const double h = 2*86400;
const double dt = 86400;
//framerate
const int framerate = 30;
const double sixth = .1666666666666666666666666666666666667;
//defines parameters of planet
class Planet {
	public:
		double posX,posY,velX,velY,r,a,ax,ay;
		void euler_step();
		void runge_step();
		void rk2_step();
};
//defines Euler Method to update position
void Planet::euler_step(){
	r = pow(pow((posX-SunX),2) + pow((posY-SunY),2),.5);
	a = (G*Ms)/(r*r);
	ax = a * (SunX - posX)/r;
	ay = a * (SunY - posY)/r;
	velX += (ax*h);
	velY += (ay*h);
	//cout << "vel " <<velX << " " << velY <<endl;
	posX += (velX*h);
	posY += (velY*h);
	cout << posX << " " << posY <<endl;

	
}
double acc(double posX,double posY,double planX, double planY){
	double r = pow((posX-planX),2) + pow((posY-planY),2);
	return (G*Ms)/(r);
}
double accX(double posX,double posY,double planX, double planY,double a){
	double n = (posX - planX)*pow(pow((posX-planX),2) + pow((posY-planY),2),-.5);
	return n*a;
}
double accY(double posX,double posY,double planX, double planY, double a){
	double n = (posY - planY)*pow(pow((posX-planX),2) + pow((posY-planY),2),-.5);
	return n*a;
}



double rad(double posX,double posY,double planX, double planY){
	return pow(pow((posX-planX),2) + pow((posY-planY),2),.5);
}

//vectorize params
vector<double> f1(vector<double> pos,vector<double> vel,vector<double> plan){
	return vel;
}

vector<double> f2(vector<double> pos,vector<double> vel,vector<double> plan){
	vector<double> vprime;
	vprime.push_back(-1*G*Ms*(pos[0]-plan[0])*pow((pow((pos[0]-plan[0]),2) + pow((pos[1]-plan[1]),2)),-1.5));
	vprime.push_back(-1*G*Ms*(pos[1]-plan[1])*pow((pow((pos[0]-plan[0]),2) + pow((pos[1]-plan[1]),2)),-1.5));
	return vprime;
}

vector<double> va(vector<double> a, vector<double> b){
	vector<double> z;
	z.push_back(a[0]+b[0]);z.push_back(a[1]+b[1]);
	return z;
}
vector<double> va4(vector<double> a, vector<double> b,vector<double> c, vector<double> d){
	vector<double> z;
	z.push_back(a[0]+b[0]+c[0]+d[0]);z.push_back(a[1]+b[1]+c[1]+d[1]);
	return z;
}

vector<double> vm(vector<double> a, double b){
	vector<double> z;
	z.push_back(a[0]*b);z.push_back(a[1]*b);
	return z;
}

void Planet::runge_step(){
	//r = pow(pow((posX-SunX),2) + pow((posY-SunY),2),.5);
	//kv1
	r = pow(pow((posX-SunX),2) + pow((posY-SunY),2),.5);
	vector<double> pos; pos.push_back(posX); pos.push_back(posY); //pos x1
	vector<double> vel; vel.push_back(velX); vel.push_back(velY); //vel x2
	vector<double> plan; plan.push_back(SunX); plan.push_back(SunY);
	vector<double> k11,k21,k12,k22,k13,k23,k14,k24;
	k11 = vm(f1(pos,vel,plan),dt);
	k21 = vm(f2(pos,vel,plan),dt);
	k12 = vm(f1(va(pos,vm(k11,.5)),va(vel,vm(k21,.5)),plan),dt);
	k22 = vm(f2(va(pos,vm(k11,.5)),va(vel,vm(k21,.5)),plan),dt);
	k13 = vm(f1(va(pos,vm(k12,.5)),va(vel,vm(k22,.5)),plan),dt);
	k23 = vm(f2(va(pos,vm(k12,.5)),va(vel,vm(k22,.5)),plan),dt);
	k14 = vm(f1(va(pos,k13),va(vel,k23),plan),dt);
	k24 = vm(f2(va(pos,k13),va(vel,k23),plan),dt);	
	pos=va(pos,vm(va4(k11,vm(k12,2),vm(k13,2),k14),sixth));
	vel=va(vel,vm(va4(k21,vm(k22,2),vm(k23,2),k24),sixth));
	
	posX = pos[0];posY = pos[1];
	velX = vel[0];velY = vel[1];
	cout << posX << " " << posY <<endl;
	
}

//average point x-----y

void Planet::rk2_step(){
	double r1 = pow(pow((posX-SunX),2) + pow((posY-SunY),2),.5);
	double a1 = (G*Ms)/(r1*r1);
	double ax1 = a1 * (SunX - posX)/r1;
	double ay1 = a1 * (SunY - posY)/r1;
	double velX1 = velX + (ax1*h);
	double velY1 = velY + (ay1*h);
	//cout << "vel " <<velX << " " << velY <<endl;
	double posX1 = posX + (velX1*h);
	double posY1 = posY + (velY1*h);
	
	double r2 = pow(pow((posX1-SunX),2) + pow((posY1-SunY),2),.5);
	double a2 = (G*Ms)/(r2*r2);
	double ax2 = a2 * (SunX - posX1)/r2;
	double ay2 = a2 * (SunY - posY1)/r2;
	double velX2 = velX + (ax2*h);
	double velY2 = velY + (ay2*h);
	//cout << "vel " <<velX << " " << velY <<endl;
	double posX2 = posX + (velX2*h);
	double posY2 = posY + (velY2*h);
	
	velX = (velX1 + velX2)*.5;
	velY = (velY1 + velY2)*.5;
	//cout << "vel " <<velX << " " << velY <<endl;
	posX += velX*h;
	posY += velY*h;
	
	
}



//assuming Msun >> Msat
int main(){
	freopen("output.txt", "w", stdout);
	//define Earth
	Planet E;
	E.posX = WinX/2 + Xi;
	E.posY = WinY/2;
	E.velX = 0;
	E.velY = Vi;
	
	//SFML Window
	 // create the window
    sf::RenderWindow window(sf::VideoMode(WinX, WinY), "EEgrav");
	window.setFramerateLimit(framerate);
    // run the program as long as the window is open
	int time = 0;
	int rad = 50;
	//define sun visual
	sf::CircleShape Sun(2*rad);
	Sun.setFillColor(sf::Color(252, 212, 64));
	Sun.setPosition(sf::Vector2f(SunX,SunY));
	//define earth visual
	//sf::CircleShape Earth(rad/5);
	sf::CircleShape Earth(rad/5);
	Earth.setFillColor(sf::Color(0, 119, 190));
	Earth.setPosition(sf::Vector2f(SunX + Xi,SunY));
	//display text
	sf::Font font;
	
	//initialize distance text
	sf::Text text_distance;
	text_distance.setFont(font); 
	text_distance.setString("Distance: " + to_string(E.r));
	text_distance.setCharacterSize(64);
	text_distance.setFillColor(sf::Color::White);
	text_distance.setStyle(sf::Text::Bold);
	//initialize perihelion text
	sf::Text text_perihelion;
	text_perihelion.setFont(font); 
	text_perihelion.setString("Perihelion: " + to_string(E.r));
	text_perihelion.setCharacterSize(64);
	text_perihelion.setFillColor(sf::Color::White);
	text_perihelion.setStyle(sf::Text::Bold);
	text_perihelion.setPosition(sf::Vector2f(0,80));
	//initialize aphelion text
	sf::Text text_aphelion;
	text_aphelion.setFont(font); 
	text_aphelion.setString("Aphelion: " + to_string(E.r));
	text_aphelion.setCharacterSize(64);
	text_aphelion.setFillColor(sf::Color::White);
	text_aphelion.setStyle(sf::Text::Bold);
	text_aphelion.setPosition(sf::Vector2f(0,160));
	//initialize year count text
	sf::Text text_year;
	text_year.setFont(font); 
	text_year.setString("Aphelion: " + to_string(E.r));
	text_year.setCharacterSize(64);
	text_year.setFillColor(sf::Color::White);
	text_year.setStyle(sf::Text::Bold);
	text_year.setPosition(sf::Vector2f(0,240));
	
	
	
	
	double aphelion = 0;
	double perihelion = pow(10,100);
	
if (!font.loadFromFile("arial.ttf"))
{
	cout << "font loading error"; 
}

    while (window.isOpen())
    {
		if(!trace){window.clear(sf::Color::Black);}		
	
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }
		//cout << E.posX << " "<< E.posY <<endl;

		
		if(runge == "rk1"){E.euler_step();}
		if(runge == "rk2"){E.rk2_step();}
		if(runge == "rk4"){E.runge_step();}

		time++;
		 
		int pixX = (int) E.posX;
		int pixY = (int) E.posY;
		if (abs(pixY - SunY)<1 && (pixX - SunX)> 0){
			text_year.setString("Year Length: "+ to_string(time*h/31536000));

		}
		Earth.setPosition(sf::Vector2f(pixX,pixY));
		//get aphelion and perihelion
		
		if(E.r > aphelion){
			aphelion = E.r;
		}
		if(E.r < perihelion){
			perihelion = E.r;
		}
		//text_year.setString("Year Length: "+ to_string(time*h/31536000));
		text_distance.setString("Distance: " + to_string(E.r * 300000) + " km (" + to_string(E.r * 300000/149598000) + " AU)");
		text_aphelion.setString("Aphelion: " + to_string(aphelion * 300000) + " km ("+ to_string(aphelion * 300000/149598000) + " AU)");
		text_perihelion.setString("Perihelion: " + to_string(perihelion * 300000) + " km (" + to_string(perihelion * 300000/149598000) + " AU)");

		
		
        window.draw(Sun);
		window.draw(Earth);	
		window.draw(text_distance);		
		window.draw(text_aphelion);		
		window.draw(text_perihelion);		
		window.draw(text_year);		

        window.display();
    }

	return 0;
}