import java.awt.BasicStroke;
import java.awt.Color;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import javax.sound.sampled.AudioSystem;
import javax.sound.sampled.Clip;

import org.opensourcephysics.controls.AbstractSimulation;
import org.opensourcephysics.controls.SimulationControl;
import org.opensourcephysics.display.Trail;
import org.opensourcephysics.frames.DisplayFrame;
import org.opensourcephysics.frames.PlotFrame;

/**
 * The BasicSpring controls the creation of a spring, oscillation, and management of the masses. 
 * It creates graphs to plot the x and y positions of masses in the spring and the 2D frame containing 
 * the spring's movement. 
 * 
 * The class uses the mouse position to control the spring, finds the proper oscillation given an amplitude 
 * and a frequency, and can perform the slinky drop. 
 * 
 * Bonuses: The graphs visually demonstrate the spring's position using a variety of colors. The mouse can be 
 * used to drag the spring, wind can be added to the simulation, it can oscillate at various frequencies, and 
 * the simulation can create an inclined plane with a large mass to demonstrate the impact of friction on the 
 * spring. It can also play a tone based on the movement of the spring. 
 * 
 * A boing sound is played when the spring changes direction, and the masses are connected by a trail. 
 * During the oscillation, the user can select which spring will be oscillated. 
 * 
 * During oscillation, the user can select a note (or amplitude/frequency) to oscillate at, or it can go through 
 * a scale by oscillating at the frequencies associated with each note. 
 * 
 * @method sumForces
 * 	sumForces takes an ArrayList of forces and calculates the sum of them (scalar sum). To find the angle, 
 * 	the calcAng method can be used.  
 * 
 * @method dist  
 * 	Calculates the distance between two particles (masses in the spring). 
 * 
 * @method calcAng
 * 	Finds the angle between two masses. The first argument is used as the center, and the angle is found with the second argument. 
 * 
 * @method play 
 * 	Plays a given .wav file. 
 * 
 * @author Andrew M. 
 */

public class BasicSpring extends AbstractSimulation {
	static DisplayFrame frame = new DisplayFrame("x", "y", "Frame");  //x-y, 2D pos frame
	static PlotFrame plotX = new PlotFrame("t", "X", "Plot X"); //plot x positions
	static PlotFrame plotY = new PlotFrame("t", "Y", "Plot Y"); //plot y positions
	ArrayList<Mass> masses = new ArrayList<Mass>();  //ArrayList of masses that comprise spring
	Trail string = new Trail(); //connect spring masses 
	int springNum = 0; //number of masses in spring
	double restLength = 0; 
	double initX = 0; 
	double mass; //mass (set by user in control) 
	double timeStep = 0.01; //default time step 
	double k; //spring constant
	double largeMass; //large mass on end for bungee
	double amplitude; //height used to drive mass making wave 
	double lowestPoint = 0; //finds lowest point of mass 
	double initLow; 
	double osc_theta; //theta to oscillate at 
	double count = 0; //count of steps (for boing sound)
	double initLength; //initial length 
	double theta_start; //starting theta 
	double magnitude, period; //for oscillation 
	double drive_mass = 0; //mass to oscillate 
	public static double g = 9.803; //gravity constant
	public static double init_deltaT; //time step
	public double alpha = 0.0000; //air resistance constant 
	public String osc_note;  
	public boolean playNote = false; 
	HashMap<Integer, Color> map = new HashMap<Integer, Color>(); //used to set the color of each mass

	boolean hasGravity = true; //whether there is gravity
	boolean fix_end; //fixed end
	boolean double_pend; 
	boolean scale; 
	double stepsPerDisplay; 
	double counter = 0; 

	protected void doStep(){
		super.setStepsPerDisplay(3); //increases speed of simulation 
		super.delayTime = 0; //end delay time 

		for (int ii = 0; ii < masses.size(); ii++) { //step through arraylist of masses
			masses.get(ii).deltaT = timeStep; //set the time step (important if scaling time) 
			masses.get(ii).forces.clear(); //clear old forces
			masses.get(ii).forces.removeAll(masses.get(ii).forces); 
			//add force of gravity (weight)
			masses.get(ii).forces.add(new Force(masses.get(ii).mass*g, -90, "mg of spring " + ii)); 

			//needed to calculate upward kx and downward kx
			addForces(ii); 
			//sum forces
			double xSum = sumForces(masses.get(ii).forces, 0); 
			if(Math.abs(xSum) < 0.005) //gets rid of imperfections in rounding with trig 
				xSum = 0; 
			double ySum = sumForces(masses.get(ii).forces, 1);
			if(Math.abs(ySum) < 0.005) //gets rid of imperfections in rounding with trig 
				ySum = 0; 
			//divide by mass to get acceleration (Newton's 2nd Law)

			//finds if mouse is within important coordinates
			masses.get(ii).acc_x = xSum/masses.get(ii).mass; //set accelerations 
			masses.get(ii).acc_y = ySum/masses.get(ii).mass;
			if(ii == 0){
				masses.get(ii).acc_x = 0; //fix first particle 
				masses.get(ii).acc_y = 0; 
			}

		}
		for (int ii = 0; ii < masses.size(); ii++) {
			masses.get(ii).Step(frame, true); //must step each simulated mass (to refresh frame) 
		}

		string.clear(); 
		for (int ii = 0; ii < masses.size(); ii++) {
			string.addPoint(masses.get(ii).getX_pos(), masses.get(ii).getY_pos()); //connect the masses
		}
		string.setStroke(new BasicStroke(2));
		frame.addDrawable(string); 
		int index = masses.size()-1; 
		plotY.append(index, masses.get(index).time, masses.get(index).getY());
		plotX.append(index, masses.get(index).time, masses.get(index).getX());
	}

	public void initialize(){
		//clear frame 
		frame.clearData();
		frame.setSize(500, 500); 
		frame.clearDrawables();
		masses.removeAll(masses); 
		masses.clear();

		plotY.setLocation(400, 0); //set location of graphs
		plotX.setLocation(800, 0); 
		//all initial values in control 
		springNum = (int) control.getDouble("Number of masses"); //set number of springs
		restLength = control.getDouble("Rest Length")/(springNum-1); //fix rest length
		initX = control.getDouble("Fixed End Location"); //set fixed x coordinate 
		fix_end = control.getBoolean("Fixed End"); //fix the end (boolean) 
		mass = control.getDouble("Mass of Spring"); //set mass 
		timeStep = control.getDouble("Timestep"); //set delta t 
		init_deltaT = timeStep; 
		k = control.getDouble("Spring Constant"); //spring constant  
		largeMass = control.getDouble("Large Mass"); //mass on end 

		g = control.getDouble("G"); //value of gravity 
		initLength = control.getDouble("Initial Length"); //initial length of spring  
		theta_start = control.getDouble("Start Theta"); //theta to hang spring at 
		//set radius of particles, add to the spring array list
		for (int ii = 0; ii < springNum; ii++) {
			masses.add(new Mass(0/*time*/)); 
			masses.get(ii).pixRadius = 3; 
			masses.get(ii).mass = mass; 
		}

		for (int ii = 0; ii < springNum; ii++) {
			if(ii == 0) { //initialize particles 
				masses.get(ii).init(initX, 0, 0, 0, 0, 0, mass, timeStep, alpha);
			}
			else {
				double xComp = initX + Math.cos(Math.toRadians(theta_start))*initLength*((ii-1)/((double)springNum-1)); //x component of init position 
				double yComp = Math.sin(Math.toRadians(theta_start))*initLength*((ii-1)/((double)springNum-1)); //y component of init position 
				masses.get(ii).init(xComp, yComp, 0, 0, 0, 0, mass, timeStep, alpha); //initialize them in initial position 
			}
			masses.get(ii).useRiemann = false; 
			frame.addDrawable(masses.get(ii));
			masses.get(ii).orig_x = masses.get(ii).x_pos; 
		}

		if(largeMass != mass)
			masses.get(masses.size()-1).setRadius(5); //just make it easier to see!
		masses.get(masses.size()-1).color = Color.GRAY.brighter(); 
		frame.setVisible(true); 

		//fill the hash map
		map.put(0, Color.red); 
		map.put(1, Color.green); 
		map.put(2, Color.blue); 
		map.put(3, Color.yellow.darker()); 
		map.put(4, Color.cyan.darker());
		map.put(5, Color.pink.brighter());
		map.put(6, Color.blue);
		map.put(7, Color.orange);
		map.put(8, Color.green.brighter());
		map.put(9, Color.blue.brighter());
		map.put(10, Color.pink.darker());
		for (int ii = 0; ii < masses.size(); ii++) {
			masses.get(ii).color = map.get(ii%map.keySet().toArray().length); //sets color of each mass (coordinate with graph) 
		}
		frame.clearDrawables(); //or clear, set gravity 
		hasGravity = true; 
		initLow = masses.get(springNum-1).y_pos; 
		masses.get(springNum-1).mass = largeMass;
	}

	public void reset(){
		//reset all conditions, location of graphs 
		plotY.setLocation(400, 0);
		plotX.setLocation(800, 0);
		control.setValue("Fixed End", true); //fixed end is default 
		control.setValue("Mouse Drag", false); 
		control.setValue("Oscillator", true); //oscillator is default 
		//all initial conditions
		control.setValue("Timestep", 0.03); 
		control.setValue("Fixed End Location", 0); 
		control.setValue("Rest Length", 0); 
		control.setValue("Initial Length", 2); 
		control.setValue("Number of masses", 30); 
		control.setValue("Mass of Spring", 1);
		control.setValue("Spring Constant", 500); 
		control.setValue("Large Mass", 6); 
		control.setValue("G", 9.803); 
		control.setValue("Start Theta", -90);

		frame.clearData();
		frame.clearDrawables();
		masses.removeAll(masses); 
		masses.clear();
	}

	public static void main(String[] args) {
		@SuppressWarnings("unused")
		SimulationControl control = SimulationControl.createApp(new BasicSpring()); //creates simulation implementing this class
	}

	/**
	 * Sums an ArrayList of forces. This method is used to calculate 
	 * acceleration on a particular mass: each mass has an ArrayList of 
	 * forces such as its weight and spring forces that must be summed 
	 * to determine acceleration (Newton's 2nd Law). 
	 * 
	 * @param f
	 * 	Arraylist of forces. 
	 * @param dir
	 * 	X or Y direction (changes whether to use x or y components of forces). 
	 * @return
	 * 	The net force in a direction (can be +/-). 
	 */
	static double sumForces(ArrayList<Force> f, double dir){
		double sum = 0; 
		if(dir == 0){
			for (int ii = 0; ii < f.size(); ii++) {
				sum += f.get(ii).forceX(); //get x components
			}
		}
		else {
			for (int ii = 0; ii < f.size(); ii++) {
				sum += f.get(ii).forceY(); //get y components
			}
		}
		return sum; 
	}

	/**
	 * Calculate the distance between two masses using distance formula. 
	 * 
	 * @param a
	 * 	Center mass. 
	 * @param b
	 * 	Target mass. 
	 * @return
	 * 	Distance between the two. 
	 */
	double dist(Mass a, Mass b){
		double dist = Math.sqrt(Math.pow(a.getX_pos() - b.getX_pos(), 2) + Math.pow(a.getY_pos() - b.getY_pos(), 2)); 
		return dist; //distance formula
	}

	/**
	 * Calculate the angle between two masses in degrees. 
	 * 
	 * @param a
	 * 	Center mass. 
	 * @param b
	 * 	Target mass. 
	 * @return
	 * 	Angle (in degrees) between the two masses. 
	 */
	public static double calcAng(Mass a, Mass b) {
		double xDiff = b.x_pos - a.x_pos; //find delta x
		double yDiff = b.y_pos - a.y_pos; //find delta y
		double angle = Math.toDegrees(Math.atan2(yDiff, xDiff));
		//angle = Math.round(angle*100)/100; 
		return angle; //inverse tangent of y/x
	}

	/**
	 * Play a sound! 
	 * 
	 * @param filename
	 * 	The name of the sound file. 
	 * 
	 * Credit: http://stackoverflow.com/questions/2416935/how-to-play-wav-files-with-java 
	 */
	public static void play(String filename){
		try {
			Clip clip = AudioSystem.getClip(); //make a new clip
			clip.open(AudioSystem.getAudioInputStream(new File(filename))); 
			clip.start(); //play clip 
		}
		catch (Exception exc){
			exc.printStackTrace(System.out);
		}
	}

	/**
	 * Adds all forces on a particular mass - wind, upward kx, downward kx, and the mass. 
	 * It calculates the correct angle and distance for the force. 
	 * 
	 * @param ii
	 * 	The number of the mass used to sum forces. 
	 */
	public void addForces(int ii){
		double upAng, downAng; //angle with mass above, mass below
		double upMag = 0, downMag = 0; //magnitude of force above, below

		//add downward kx
		if(ii != springNum-1){ //the last spring does not have a downward kx
			double dist = dist(masses.get(ii), masses.get(ii+1)); //distance for spring force

			downMag = k*(dist-restLength); //calculate kx
			downAng = calcAng(masses.get(ii), masses.get(ii+1)); //calculate angle between two masses 
			masses.get(ii).forces.add(new Force(downMag, downAng, "DOWN kx between spring " + ii + " and " + (ii+1))); //add as force
		}
		//add upward kx
		if(ii != 0){ //the first spring does not have an upward kx 
			double dist = dist(masses.get(ii), masses.get(ii-1)); //distance for spring force

			upMag = k*(dist-restLength); //calculate kx
			upAng = calcAng(masses.get(ii), masses.get(ii-1)); //calculate angle between two masses
			masses.get(ii).forces.add(new Force(upMag, upAng, "UP kx between spring " + ii + " and " + (ii-1))); //add as force
		}
	}
}
