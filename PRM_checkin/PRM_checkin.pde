import java.util.TreeMap; //<>//
import java.util.Iterator;

class Object{
  PVector pos;
  Object(PVector l) {
    pos = l.copy();
    }
  void draw_obj(){
    noStroke();
    fill(0,255,0);
    translate(pos.x,0,pos.y);
    drawCylinder(10, 10, cylinder_height, 64); 
    translate(-pos.x,0,-pos.y);
    fill(255);
  }
}

boolean is_collsion(PVector point1, PVector point2, PVector origin, float radius){
  float x1 = point1.x - origin.x;
  float y1 = point1.y - origin.y;
  float x2 = point2.x - origin.x;
  float y2 = point2.y - origin.y;
  float dx = x2-x1;
  float dy = y2-y1;
  float dr = sqrt(dx*dx+dy*dy);
  float D = x1*y2-x2*y1;
  float delta = radius*radius*dr*dr - D*D;
  if(delta<=0){
    return false;
  }else{
    float y0 = (-D*dx+dy*sqrt(delta))/(dr*dr);
    float ymul = (y0-y1)*(y0-y2);
    if (ymul>0){
      return false;
    }else{
      return true;
    }
  }
}

PVector get_closest_point_on_circle(PVector point, PVector origin,float radius){
  float Px = point.x - origin.x;
  float Py = point.y - origin.y;
  float tmp = sqrt(Px*Px*radius*radius/(Px*Px+Py*Py));
  float sgn = 1;
  if (Px < 0) {
    sgn = -1;
  }
  float tmp_x = sgn*tmp;
  float tmp_y = (Py/Px)*tmp_x;
  float final_x = tmp_x + origin.x;
  float final_y = tmp_y + origin.y;
  return new PVector(final_x,final_y);
}

void drawCylinder(float topRadius, float bottomRadius, float tall, int sides) {
  float angle = 0;
  float angleIncrement = TWO_PI / sides;
  beginShape(QUAD_STRIP);
  for (int i = 0; i < sides + 1; ++i) {
    vertex(topRadius*cos(angle), 0, topRadius*sin(angle));
    vertex(bottomRadius*cos(angle), tall, bottomRadius*sin(angle));
    angle += angleIncrement;
  }
  endShape();
  if (topRadius != 0) {
    angle = 0;
    beginShape(TRIANGLE_FAN);
    vertex(0, 0, 0);
    for (int i = 0; i < sides + 1; i++) {
      vertex(topRadius * cos(angle), 0, topRadius * sin(angle));
      angle += angleIncrement;
    }
    endShape();
  }
}

void draw_floor(){
  fill(0,0,255);
  rotateX(PI/2);
  translate(0,cylinder_height+1,0);
  rect(-100,-100, 200, 200);
  rotateX(-PI/2);
  translate(0,-cylinder_height-1,0);
  fill(0);
}

void draw_obstacle(){
  fill(255);
  drawCylinder(20, 20, cylinder_height, 64);
  fill(0);
}

Object gamechar = new Object(new PVector(-90,-90));

//ArrayList<Vertex> milestone = new ArrayList<Vertex>();
//ArrayList<Edge> edges = new ArrayList<Edge>();

Vertex startpoint;
Vertex endpoint;
Graph graph=new Graph();
ArrayList<Vertex> vertexs;
ArrayList<Edge> edges;
ArrayList<Vertex> path;

Camera camera;
float cylinder_height = 50;
int start_num;

void add_vertex(){
  graph.addVertex(startpoint);
  graph.addVertex(endpoint);
  //milestone.add(startpoint);
  //milestone.add(endpoint);
  int num_midpoint = 2;
  //randomSeed(5611);
  while(num_midpoint<100){
    PVector new_l = new PVector(random(-90,90),random(-90,90));
    float dis_to_obs1 = PVector.dist(new_l, new PVector(0,0));
    if (dis_to_obs1>=30){
      graph.addVertex(new Vertex(new_l,num_midpoint));
      //milestone.add(new Vertex(new_l));
      num_midpoint = num_midpoint +1;
    }else{
      //here we add 5 to 30 so that the closest point not collsion
      PVector closest_l = get_closest_point_on_circle(new_l,new PVector(0,0),30+5);
      graph.addVertex(new Vertex(closest_l,num_midpoint));
      num_midpoint = num_midpoint +1;
    }
  }
}

void add_edge(){
  for (int i = 0; i < vertexs.size(); i++) {
    Vertex point1 = vertexs.get(i);
    ArrayList<Float> dist = new ArrayList<Float>();
    for (int j = 0; j < vertexs.size(); j++) {
      Vertex point2 = vertexs.get(j);
      dist.add(-point2.pos_2D.dist(point1.pos_2D)); 
    }
    float[] floatArray = new float[dist.size()];
    int p = 0;
    for (Float f : dist) {
      floatArray[p++] = (f != null ? f : Float.NaN); // Or whatever default you want.
    }
    int[] rst  = indexesOfTopElements(floatArray,6);
    for (int k = 1; k < 6; k++) {
      Vertex p1 = point1;
      Vertex p2 = vertexs.get(rst[k]);
      if(!is_collsion(p1.pos_2D,p2.pos_2D,new PVector(0,0),30)){
        Edge new_edge =new Edge(p1,p2);
        if (new_edge.dist<50){
          graph.addEdge(new_edge);
          graph.addEdge(new Edge(p2,p1));
        }
      }
    }
  }
}

boolean move_from_to(Vertex vertex1, Vertex vertex2,float dt){
  boolean is_arrive = false;
  PVector pos1 = vertex1.pos_2D;
  PVector pos2 = vertex2.pos_2D;
  gamechar.pos = pos1;
  //use mult to contorl speed
  PVector spd = pos2.copy().sub(pos1).normalize().mult(20);
  float dist_to_dest = PVector.dist(gamechar.pos, pos2);
  //println(dist_to_dest);
  if(dist_to_dest>0.5){
    gamechar.pos = gamechar.pos.add(spd.mult(dt));
    gamechar.draw_obj();
    return is_arrive;
  }else{
    gamechar.draw_obj();
    is_arrive = true;
    return is_arrive;
  }
}

void setup(){
  size(1000, 750, P3D);
  surface.setTitle("Motion Planning");
  camera = new Camera();
  noStroke();   
  noFill();
  startpoint = new Vertex(new PVector(-90,-90),0);
  endpoint = new Vertex(new PVector(90,90),1);
  add_vertex();
  vertexs=graph.vertexs;
  add_edge();
  edges=graph.edges;
  println(edges.size()/2);
  path = BFS(graph);
  if(path.size() == 0){
    gamechar = new Object(new PVector(-90,-90));
    graph=new Graph();
    vertexs.get(0).visited = false;
    vertexs.get(1).visited = false;
    setup();
  }
  println(path.size());
  start_num = path.size()-1;

}

void draw_current_pos(){
  gamechar.draw_obj();
}

void draw(){
  float startFrame = millis();
  background(0);
  lights();
  camera.Update(1.0/frameRate);
  draw_floor();
  draw_obstacle();
  for (int i = 0; i < vertexs.size(); i++) {
    Vertex midpoint = vertexs.get(i);
    midpoint.draw_vert();
  }
  for (int i = 0; i < edges.size(); i++) {
    Edge midedge = edges.get(i);
    midedge.draw_edge();
  }
  float endPhysics = millis();
  
  float dt = 1.0/frameRate;

  if (start_num>0){
    Vertex curr = path.get(start_num);
    Vertex next = path.get(start_num-1);
    boolean arrive_cond = move_from_to(curr, next, dt);
    if (arrive_cond){
      start_num = start_num-1;
    }
  }else{
    Vertex curr = path.get(0);
    Vertex next = endpoint;
    boolean arrive_cond = move_from_to(curr, next, dt);
    if (arrive_cond){
      draw_current_pos();
    }
  }
  
  //for (int i = path.size()-1; i >= 0; i--) {
  //  Vertex curr = path.get(i);
  //  //println(i,curr.index);
  //}
  
  
  
  float endFrame = millis();
  String runtimeReport = "Frame: "+str(endFrame-startFrame)+"ms,"+
        " Physics: "+ str(endPhysics-endFrame)+"ms,"+
        " FPS: "+ str(round(frameRate)) +"\n";
  surface.setTitle("PRM_check"+ "  -  " +runtimeReport);
}

void keyPressed()
{
  camera.HandleKeyPressed();
  if (key =='r'){
    gamechar = new Object(new PVector(-90,-90));
    graph=new Graph();
    vertexs.get(0).visited = false;
    vertexs.get(1).visited = false;
    setup();
  }
  if (key =='t'){
    print(camera.position);
  }
}

void keyReleased()
{
  camera.HandleKeyReleased();
}

void mouseClicked(MouseEvent evt) {
  if (evt.getCount() == 2) doubleClicked();
}

void doubleClicked() {
  
}

  
class Camera
{
  Camera()
  {
    position      = new PVector(-400,-300,-10); // initial position
    theta         = -1.6; // rotation around Y axis. Starts with forward direction as ( 0, 0, -1 )
    phi           = -0.5; // rotation around X axis. Starts with up direction as ( 0, 1, 0 )
    moveSpeed     = 200;
    turnSpeed     = 1.57; // radians/sec
    
    // dont need to change these
    negativeMovement = new PVector( 0, 0, 0 );
    positiveMovement = new PVector( 0, 0, 0 );
    negativeTurn     = new PVector( 0, 0 ); // .x for theta, .y for phi
    positiveTurn     = new PVector( 0, 0 );
    fovy             = PI / 4;
    aspectRatio      = width / (float) height;
    nearPlane        = 0.1;
    farPlane         = 10000;
  }
  
  void Update( float dt )
  {
    theta += turnSpeed * (negativeTurn.x + positiveTurn.x) * dt;
    
    // cap the rotation about the X axis to be less than 90 degrees to avoid gimble lock
    float maxAngleInRadians = 85 * PI / 180;
    phi = min( maxAngleInRadians, max( -maxAngleInRadians, phi + turnSpeed * ( negativeTurn.y + positiveTurn.y ) * dt ) );
    
    // re-orienting the angles to match the wikipedia formulas: https://en.wikipedia.org/wiki/Spherical_coordinate_system
    // except that their theta and phi are named opposite
    float t = theta + PI / 2;
    float p = phi + PI / 2;
    PVector forwardDir = new PVector( sin( p ) * cos( t ),   cos( p ),   -sin( p ) * sin ( t ) );
    PVector upDir      = new PVector( sin( phi ) * cos( t ), cos( phi ), -sin( t ) * sin( phi ) );
    PVector rightDir   = new PVector( cos( theta ), 0, -sin( theta ) );
    PVector velocity   = new PVector( negativeMovement.x + positiveMovement.x, negativeMovement.y + positiveMovement.y, negativeMovement.z + positiveMovement.z );
    position.add( PVector.mult( forwardDir, moveSpeed * velocity.z * dt ) );
    position.add( PVector.mult( upDir,      moveSpeed * velocity.y * dt ) );
    position.add( PVector.mult( rightDir,   moveSpeed * velocity.x * dt ) );
    
    aspectRatio = width / (float) height;
    perspective( fovy, aspectRatio, nearPlane, farPlane );
    camera( position.x, position.y, position.z,
            position.x + forwardDir.x, position.y + forwardDir.y, position.z + forwardDir.z,
            upDir.x, upDir.y, upDir.z );
  }
  
  // only need to change if you want difrent keys for the controls
  void HandleKeyPressed()
  {
    if ( key == 'w' ) positiveMovement.z = 1;
    if ( key == 's' ) negativeMovement.z = -1;
    if ( key == 'a' ) negativeMovement.x = -1;
    if ( key == 'd' ) positiveMovement.x = 1;
    if ( key == 'q' ) positiveMovement.y = 1;
    if ( key == 'e' ) negativeMovement.y = -1;
    
    if ( keyCode == LEFT )  negativeTurn.x = 1;
    if ( keyCode == RIGHT ) positiveTurn.x = -1;
    if ( keyCode == UP )    positiveTurn.y = 1;
    if ( keyCode == DOWN )  negativeTurn.y = -1;
  }
  
  // only need to change if you want difrent keys for the controls
  void HandleKeyReleased()
  {
    if ( key == 'w' ) positiveMovement.z = 0;
    if ( key == 'q' ) positiveMovement.y = 0;
    if ( key == 'd' ) positiveMovement.x = 0;
    if ( key == 'a' ) negativeMovement.x = 0;
    if ( key == 's' ) negativeMovement.z = 0;
    if ( key == 'e' ) negativeMovement.y = 0;
    
    if ( keyCode == LEFT  ) negativeTurn.x = 0;
    if ( keyCode == RIGHT ) positiveTurn.x = 0;
    if ( keyCode == UP    ) positiveTurn.y = 0;
    if ( keyCode == DOWN  ) negativeTurn.y = 0;
  }
  
  // only necessary to change if you want different start position, orientation, or speeds
  PVector position;
  float theta;
  float phi;
  float moveSpeed;
  float turnSpeed;
  
  // probably don't need / want to change any of the below variables
  float fovy;
  float aspectRatio;
  float nearPlane;
  float farPlane;  
  PVector negativeMovement;
  PVector positiveMovement;
  PVector negativeTurn;
  PVector positiveTurn;
};
