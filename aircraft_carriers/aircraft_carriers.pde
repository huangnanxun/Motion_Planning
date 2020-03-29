import java.util.TreeMap; //<>//
import java.util.Iterator;

class Object{
  PVector pos;
  Object(PVector l) {
    pos = l.copy();
    }
  void draw_obj(int obj_num){
    noStroke();
    fill(0,255,0);
    translate(pos.x,cylinder_height,pos.y);
    //drawCylinder(10, 10, cylinder_height, 64);
    if (obj_num == 1){
    rotateX(-PI);
    rotateY(-PI);
    ship.setFill(color(255,0,0));
    shape(ship);
    rotateY(PI);
    rotateX(PI); 
    }
    if (obj_num == 2){
    rotateX(-PI);
    fill(0,255,0);
    ship.setFill(color(0,255,0));
    shape(ship);
    rotateX(PI);
    }
    translate(-pos.x,-cylinder_height,-pos.y);
    fill(255);
  }
}

boolean is_collsion_circle(PVector point1, PVector point2, PVector origin, float radius){
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

boolean is_collsion_rect(PVector point1, PVector point2, float rx, float ry, float rw, float rh) {
  float x1 = point1.x;
  float y1 = point1.y;
  float x2 = point2.x;
  float y2 = point2.y;
  boolean left =   lineLine(x1,y1,x2,y2, rx,ry,rx, ry+rh);
  boolean right =  lineLine(x1,y1,x2,y2, rx+rw,ry, rx+rw,ry+rh);
  boolean top =    lineLine(x1,y1,x2,y2, rx,ry, rx+rw,ry);
  boolean bottom = lineLine(x1,y1,x2,y2, rx,ry+rh, rx+rw,ry+rh);
  if (left || right || top || bottom) {
    return true;
  }
  return false;
}

float point_rec_dist(PVector point,float rx, float ry, float rw, float rh){
  float x = point.x;
  float y = point.y;
  float dx = max(abs(rx - x)-rw/2, 0);
  float dy = max(abs(ry - y)-rh/2, 0);
  return sqrt(dx*dx+dy*dy);
}


boolean lineLine(float x1, float y1, float x2, float y2, float x3, float y3, float x4, float y4) {
  float uA = ((x4-x3)*(y1-y3) - (y4-y3)*(x1-x3)) / ((y4-y3)*(x2-x1) - (x4-x3)*(y2-y1));
  float uB = ((x2-x1)*(y1-y3) - (y2-y1)*(x1-x3)) / ((y4-y3)*(x2-x1) - (x4-x3)*(y2-y1));
  if (uA >= 0 && uA <= 1 && uB >= 0 && uB <= 1) {
    float intersectionX = x1 + (uA * (x2-x1));
    float intersectionY = y1 + (uA * (y2-y1));
    fill(255,0,0);
    noStroke();
    ellipse(intersectionX, intersectionY, 20, 20);
    return true;
  }
  return false;
}


boolean check_cuurent_collsion(Object obj, PVector target_pos){
  if(is_collsion_circle(obj.pos,target_pos,new PVector(0,0),30) || is_collsion_rect(obj.pos,target_pos,-100,0,50+10,100+10)){
    return true;
  }
  return false;
}

boolean check_in_box(PVector target_pos, float box_w, float box_d, float box_dep_x, float box_dep_z){
  float lower_x = box_dep_x - 0.5*box_w - 5;
  float higher_x = box_dep_x + 0.5*box_w + 5;
  float lower_z = box_dep_z - 0.5*box_d - 5;
  float higher_z = box_dep_z + 0.5*box_d + 5;
  if (target_pos.x >lower_x && target_pos.x < higher_x && target_pos.y > lower_z && target_pos.y < higher_z){
    return true;
  }else{
    return false;
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
  rect(-200,-200, 400, 400);
  rotateX(-PI/2);
  translate(0,-cylinder_height-1,0);
  fill(0);
}

void draw_obstacle(){
  fill(255);
  drawCylinder(20, 20, cylinder_height, 64);
  fill(0);
}

void draw_obstacle2(){
  fill(255);
  translate(-100,0.5*cylinder_height,0);
  box(50, cylinder_height, 100);
  translate(100,-0.5*cylinder_height,0);
  fill(0);
}


Object gamechar = new Object(new PVector(-190,-190));
Object gamechar2 = new Object(new PVector(-190,190));
Object gamechar3 = new Object(new PVector(190,-190));

//ArrayList<Vertex> milestone = new ArrayList<Vertex>();
//ArrayList<Edge> edges = new ArrayList<Edge>();

Vertex startpoint;
Vertex endpoint;
Vertex startpoint2;
Vertex endpoint2;

Vertex obstacle2_point1;
Vertex obstacle2_point2;
Vertex obstacle2_point3;
Vertex obstacle2_point4;


Graph graph=new Graph();

ArrayList<Vertex> vertexs;
ArrayList<Edge> edges;

ArrayList<Vertex> path;
ArrayList<Vertex> path2;
ArrayList<Vertex> path3;

PShape ship;

Flock flock;

Camera camera;
float cylinder_height = 50;
int start_num;
int start_num2;
int start_num3;

void add_vertex(){
  graph.addVertex(startpoint);
  graph.addVertex(endpoint);
  graph.addVertex(startpoint2);
  graph.addVertex(endpoint2);
  
  graph.addVertex(obstacle2_point1);
  graph.addVertex(obstacle2_point2);
  graph.addVertex(obstacle2_point3);
  graph.addVertex(obstacle2_point4);
  
  randomSeed(2);
  int num_midpoint = 8;
  // 1 2 108
  
  while(num_midpoint<108){
    PVector new_l = new PVector(random(-190,190),random(-190,190));
    float dis_to_obs1 = PVector.dist(new_l, new PVector(0,0));
    if (dis_to_obs1>=30 && !check_in_box(new_l,50,100,-100,0)){
      graph.addVertex(new Vertex(new_l,num_midpoint));
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
      floatArray[p++] = (f != null ? f : Float.NaN);
    }
    int[] rst  = indexesOfTopElements(floatArray,6);
    for (int k = 1; k < 6; k++) {
      Vertex p1 = point1;
      Vertex p2 = vertexs.get(rst[k]);
      if(!is_collsion_circle(p1.pos_2D,p2.pos_2D,new PVector(0,0),30)&&!is_collsion_rect(p1.pos_2D,p2.pos_2D,-100,0,50+10,100+10)){
        Edge new_edge =new Edge(p1,p2);
        if (new_edge.dist<500){
          graph.addEdge(new_edge);
          graph.addEdge(new Edge(p2,p1));
        }
      }
    }
  }
}

boolean move_from_to(Object obj, Object obj2, Vertex vertex1, Vertex vertex2,float dt,int obj_num){
  boolean is_arrive = false;
  PVector pos1 = vertex1.pos_2D;
  PVector pos2 = vertex2.pos_2D;
  obj.pos = pos1;
  //use mult to contorl speed
  PVector spd = pos2.copy().sub(pos1).normalize().mult(40);
  flock.run();
  float dist_to_dest = PVector.dist(obj.pos, pos2);
  //println(dist_to_dest);
  float obj_dist = PVector.dist(obj.pos,obj2.pos);
  if (obj_dist<40){
    spd.add(obj2.pos.copy().sub(obj.pos).normalize().mult(-40));
  }
  if(dist_to_dest>0.5){
    obj.pos = obj.pos.add(spd.mult(dt));
    obj.draw_obj(obj_num);
    return is_arrive;
  }else{
    obj.draw_obj(obj_num);
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
  startpoint = new Vertex(new PVector(-190,-190),0);
  endpoint = new Vertex(new PVector(190,190),1);
  startpoint2 = new Vertex(new PVector(-190,190),2);
  endpoint2 = new Vertex(new PVector(190,-190),3);
  
  obstacle2_point1 = new Vertex(new PVector(-75+5,-50-5),4);
  obstacle2_point2 = new Vertex(new PVector(-75+5,50+5),5);
  obstacle2_point3 = new Vertex(new PVector(-125-5,-50-5),6);
  obstacle2_point4 = new Vertex(new PVector(-125-5,50+5),7);
  ship = loadShape("ship.obj");
  ship.scale(0.15);
  
  add_vertex();
  vertexs=graph.vertexs;
  add_edge();
  edges=graph.edges;
  path = A_Star(graph,0,1,0);
  print("path1 created");
  for(Vertex vertex:vertexs){
    vertex.visited=false;
  }
  path2 = A_Star(graph,2,3,1);
  
  print(path.size());
  if(path.size() == 0){
    gamechar = new Object(new PVector(-190,-190));
    graph=new Graph();
    vertexs.get(0).visited = false;
    vertexs.get(1).visited = false;
    setup();
  }
  print("path2 created");
  print(path2.size());
  if(path2.size() == 0){
    gamechar2 = new Object(new PVector(-190,-190));
    graph=new Graph();
    vertexs.get(2).visited = false;
    vertexs.get(3).visited = false;
    setup();
  }
  
  //for(Vertex vertex:vertexs){
  //  vertex.visited=false;
  //}
  
  //path3 = A_Star(graph,3,1,2);
  //print(path3.size());
  //if(path3.size() == 0){
  //  gamechar2 = new Object(new PVector(190,-190));
  //  graph=new Graph();
  //  vertexs.get(3).visited = false;
  //  vertexs.get(1).visited = false;
  //  setup();
  //}
  
  start_num = path.size()-1;
  start_num2 = path2.size()-1;
  //start_num3 = path3.size()-1;
  
  PVector vel1 = endpoint.pos_2D.copy().sub(startpoint.pos_2D).normalize();
  PVector vel2 = endpoint2.pos_2D.copy().sub(startpoint2.pos_2D).normalize();
  flock = new Flock();
  for (int i = 0; i < 25; i++) {
    flock.addBoid(new Boid(startpoint.pos_2D.x,startpoint.pos_2D.y,vel1,color(255,0,0)));
  }
  for (int i = 0; i < 25; i++) {
    flock.addBoid(new Boid(startpoint2.pos_2D.x,startpoint2.pos_2D.y,vel2,color(0,255,0)));
  }

}


void draw(){
  float startFrame = millis();
  background(0);
  lights();
  camera.Update(1.0/frameRate);
  draw_floor();
  draw_obstacle();
  draw_obstacle2();
  //for (int i = 0; i < vertexs.size(); i++) {
  //  Vertex midpoint = vertexs.get(i);
  //  midpoint.draw_vert();
  //}
  //for (int i = 0; i < edges.size(); i++) {
  //  Edge midedge = edges.get(i);
  //  midedge.draw_edge();
  //}
  float endPhysics = millis();
  
  float dt = 1.0/frameRate;
  
  //move(start_num, gamechar, dt, path, endpoint);
  if (start_num>1){
    Vertex curr = path.get(start_num);
    Vertex next = path.get(start_num-1);
    Vertex next_next = path.get(start_num-2);
    boolean arrive_cond = move_from_to(gamechar,gamechar2,curr, next, dt,1);
    if (arrive_cond){
      start_num = start_num-1;
    }
    if(!check_cuurent_collsion(gamechar,next_next.pos_2D)){
      start_num = start_num-1;
      int curr_index = path.get(start_num).index;
      path.set(start_num,new Vertex(gamechar.pos,curr_index));
  //    println("change to index");
  //    print(start_num);
    }
  }else{
    if(!check_cuurent_collsion(gamechar,endpoint.pos_2D)){
      Vertex curr = new Vertex(gamechar.pos,0);
      Vertex next = endpoint;
     move_from_to(gamechar,gamechar2,curr, next, dt,1);
    }else{
      Vertex curr = new Vertex(gamechar.pos,0);
      Vertex next = path.get(0);
      move_from_to(gamechar,gamechar2,curr, next, dt,1);
    }
  }
  
  if (start_num2>1){
    Vertex curr = path2.get(start_num2);
    Vertex next = path2.get(start_num2-1);
    Vertex next_next = path2.get(start_num2-2);
    boolean arrive_cond = move_from_to(gamechar2,gamechar,curr, next, dt,2);
    if (arrive_cond){
      start_num2 = start_num2-1;
    }
    if(!check_cuurent_collsion(gamechar2,next_next.pos_2D)){
      start_num2 = start_num2-1;
      int curr_index = path2.get(start_num).index;
      path2.set(start_num2,new Vertex(gamechar2.pos,curr_index));
  //    println("change to index");
  //    print(start_num);
    }
  }else{
    if(!check_cuurent_collsion(gamechar2,endpoint2.pos_2D)){
      Vertex curr = new Vertex(gamechar2.pos,0);
      Vertex next = endpoint2;
     move_from_to(gamechar2,gamechar,curr, next, dt,2);
    }else{
      Vertex curr = new Vertex(gamechar2.pos,0);
      Vertex next = path2.get(0);
      move_from_to(gamechar2,gamechar,curr, next, dt,2);
    }
  }
  
  
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
