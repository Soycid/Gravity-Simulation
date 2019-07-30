/*
TODO:
Add moons and other planets for more accuracy
Add option for relativistic model?
Figure out way to graph 
 >super trails and multiple integrators simulatanuously for visualization
 >trace EXACT solution
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
string runge = "rk1";
bool trace = false;
bool path = true;


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
const double Me = 5.972*pow(10,24);
const double Mm = 7.348*pow(10,22);


const double scale = 300000;

//Initial Position of Earh in pix (x, 0)
//const double Xi = 147099761/scale;
//Initial Velocity of Earh in pix/s
//const double Vi = pow(G*Ms*(2/Xi - 1/((1.496/3)*pow(10,3))),.5);
const double Xi = 152041481.3/scale;
const double Yi = 0/scale;
const double Vxi = -0/scale;
const double Vyi = -29.29445088/scale;
//Initial Position of Moon in pix (x, 0)
const double mXi = 360708.9606/scale;
const double mYi = -0/scale;
const double mVxi = 0/scale;
const double mVyi = 1.090492451/scale;

//Initial Velocity of Moon in pix/s
//const double mVi = pow(G*Ms*(2/Xi - 1/((1.496/3)*pow(10,3))),.5);



const double stepSpeed = 1;
//Step size in seconds
const double dt = 86400/stepSpeed;
//framerate
const int framerate = 30*stepSpeed;
const double sixth = .1666666666666666666666666666666666667;
//defines parameters of planet
class Planet {
	public:
		double posX,posY,velX,velY,r,a,ax,ay;
		std::vector<sf::Vertex> trail;
		double uposX = 0;
		double uposY = 0;
		void euler_step(double SunX, double SunY,double Ms);
		void runge_step(double SunX, double SunY);
		void rk2_step(double SunX, double SunY);
		void verlet_step(double SunX, double SunY);
		void reuler_step(double SunX, double SunY);
		void update();

		
		
		
};

void Planet::update(){
	posX += uposX;
	posY += uposY;
	uposX = 0;
	uposY = 0;
}
//defines Euler Method to update position
void Planet::euler_step(double X, double Y,double M){
	r = pow(pow((posX-X),2) + pow((posY-Y),2),.5);
	a = (G*M)/(r*r);
	ax = a * (X - posX)/r;
	ay = a * (Y - posY)/r;
	velX += (ax*dt);
	velY += (ay*dt);
	//cout << "vel " <<velX << " " << velY <<endl;
	uposX += (velX*dt);
	uposY += (velY*dt);
	//cout << posX << " " << posY <<endl;

	
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
vector<double> va3(vector<double> a, vector<double> b,vector<double> c){
	vector<double> z;
	z.push_back(a[0]+b[0]+c[0]);z.push_back(a[1]+b[1]+c[1]);
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
void Planet::rk2_step(double SunX, double SunY){
	r = pow(pow((posX-SunX),2) + pow((posY-SunY),2),.5);
	vector<double> pos; pos.push_back(posX); pos.push_back(posY); //pos x1
	vector<double> vel; vel.push_back(velX); vel.push_back(velY); //vel x2
	vector<double> plan; plan.push_back(SunX); plan.push_back(SunY);
	vector<double> k11,k21,k12,k22;
	k11 = vm(f1(pos,vel,plan),dt);
	k21 = vm(f2(pos,vel,plan),dt);
	k12 = vm(f1(va(pos,k11),va(vel,k21),plan),dt);
	k22 = vm(f2(va(pos,k11),va(vel,k21),plan),dt);

	pos=va(pos,vm(va(k11,k12),.5));
	vel=va(vel,vm(va(k21,k22),.5));

	posX = pos[0];posY = pos[1];
	velX = vel[0];velY = vel[1];
}

void Planet::runge_step(double SunX, double SunY){
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
	//cout << posX << " " << posY <<endl;
	
}

void Planet::verlet_step(double SunX, double SunY){
	r = pow(pow((posX-SunX),2) + pow((posY-SunY),2),.5);
	vector<double> pos; pos.push_back(posX); pos.push_back(posY); //pos x1
	vector<double> pos_new;
	vector<double> vel; vel.push_back(velX); vel.push_back(velY); //vel x2
	vector<double> plan; plan.push_back(SunX); plan.push_back(SunY);
	
	pos_new = va3(pos,vm(vel,dt),vm(f2(pos,vel,plan),.5*dt*dt));
	vel = va(vel,vm(va(f2(pos,vel,plan),f2(pos_new,vel,plan)),.5*dt));
	
	posX = pos_new[0]; posY = pos_new[1];
	velX = vel[0]; velY = vel[1];
}

void Planet::reuler_step(double SunX, double SunY){
	r = pow(pow((posX-SunX),2) + pow((posY-SunY),2),.5);
	vector<double> pos; pos.push_back(posX); pos.push_back(posY); //pos x1
	vector<double> vel; vel.push_back(velX); vel.push_back(velY); //vel x2
	vector<double> plan; plan.push_back(SunX); plan.push_back(SunY);
	
	vel = va(vel,vm(f2(pos,vel,plan),dt));
	pos = va(pos,vm(vel,dt));
	
	posX = pos[0]; posY = pos[1];
	velX = vel[0]; velY = vel[1];
	
}



//assuming Msun >> Msat
int main(){
	freopen("output.csv", "w", stdout);
	//define Earth
	Planet E;
	E.posX = WinX/2 + Xi;
	E.posY = WinY/2 + Yi;
	E.velX = Vxi;
	E.velY = Vyi;
	
	Planet M;
	M.posX = WinX/2 + Xi + mXi;
	M.posY = WinY/2 + Yi + mYi;
	M.velX = Vxi + mVxi;
	//M.velY = Vi + mVi;
	M.velY = Vyi + mVyi;
	M.velY = Vyi;

	
	
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
	Sun.setPosition(sf::Vector2f(SunX-2*rad,SunY-2*rad));
	//define earth visual
	sf::CircleShape Earth(rad/5);
	Earth.setFillColor(sf::Color(0, 119, 190));
	Earth.setPosition(sf::Vector2f(SunX + Xi,SunY + Yi));
	//define moon Visual
	sf::CircleShape Moon(rad/8);
	Moon.setFillColor(sf::Color(254, 252, 252));
	Moon.setPosition(sf::Vector2f(SunX + Xi + mXi,SunY + Yi + mYi));
	

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
	//initialize energy text
	sf::Text text_energy;
	text_energy.setFont(font); 
	text_energy.setString("Energy " + to_string(0));
	text_energy.setCharacterSize(64);
	text_energy.setFillColor(sf::Color::White);
	text_energy.setStyle(sf::Text::Bold);
	text_energy.setPosition(sf::Vector2f(0,320));
	
	
	
	
	cout << "time,time(sol days),Distance from Sun (km), Distance from Sun (AU),max/min"<<endl;
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

		//step code
		//if(runge == "rk1"){E.euler_step(SunX,SunY,Ms);M.euler_step(SunX,SunY,Ms);E.euler_step(M.posX,M.posY,Mm);M.euler_step(E.posX,E.posY,Me);}
		if(runge == "rk1"){E.euler_step(SunX,SunY,Ms);M.euler_step(SunX,SunY,Ms);M.euler_step(E.posX,E.posY,Me);}
		if(runge == "rk2"){E.rk2_step(SunX,SunY);}
		if(runge == "rk4"){E.runge_step(SunX,SunY);}
		if(runge == "reuler"){E.reuler_step(SunX,SunY);}
		if(runge == "verlet"){E.verlet_step(SunX,SunY);}
		E.update();
		M.update();
		//output intit
		cout << time << "," << time*dt/86400 << "," << E.r * scale <<"," << E.r * scale/149598000<< ","; 

		
		time++;
		 //energy
		 double energy = .5*(pow(E.velX,2)+pow(E.velY,2))*Me - G*Ms*Me/(E.r*E.r);
		 
		 
		
		//trail code
		if(path){
			sf::Vertex vertex;
			//changed pix to pos
			vertex.position = sf::Vector2f(E.posX, E.posY);
			vertex.color = sf::Color::Red;
			E.trail.push_back(vertex);
		}
		//change
		if (abs(E.posY - SunY)<1){
			text_year.setString("Year Length: "+ to_string(time*dt/31536000));

		}
		Earth.setPosition(sf::Vector2f(E.posX,E.posY));
		Moon.setPosition(sf::Vector2f(M.posX,M.posY));

		//get aphelion and perihelion
		
		if(E.r > aphelion){
			aphelion = E.r;
			cout << "aphelion,";
			cout << energy;
		}
		if(E.r < perihelion){
			perihelion = E.r;
			cout << "perihelion";
		}
		//text_year.setString("Year Length: "+ to_string(time*dt/31536000));
		text_distance.setString("Distance: " + to_string(E.r * scale) + " km (" + to_string(E.r * scale/149598000) + " AU)");
		text_aphelion.setString("Aphelion: " + to_string(aphelion * scale) + " km ("+ to_string(aphelion * scale/149598000) + " AU)");
		text_perihelion.setString("Perihelion: " + to_string(perihelion * scale) + " km (" + to_string(perihelion * scale/149598000) + " AU)");
		text_energy.setString("Energy: " + to_string(energy));

		
		if(path){window.draw(&E.trail[0], E.trail.size(), sf::Lines);}
        window.draw(Sun);
		window.draw(Earth);
		window.draw(Moon);	
		
		window.draw(text_distance);		
		window.draw(text_aphelion);		
		window.draw(text_perihelion);		
		window.draw(text_year);
		window.draw(text_energy);		
		
		cout << endl;
        window.display();
    }

	return 0;
}