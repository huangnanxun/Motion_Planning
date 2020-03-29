//import java.util.Comparator; //<>// //<>//
//import java.util.Map;
//import java.util.Map.Entry;
import java.util.Collections;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.Queue;

class Graph{
  ArrayList<Vertex> vertexs=new ArrayList<Vertex>();
  ArrayList<Edge> edges=new ArrayList<Edge>();
  public void addVertex(Vertex vertex) {
    vertexs.add(vertex);
  }
  public void addEdge(Edge edge) {
      edges.add(edge);
      }
}

class Vertex{
  PVector pos_2D;
  boolean visited=false;
  Vertex current_pre0;
  Vertex current_pre1;
  Vertex current_pre2;
  int index;
  Vertex(PVector l,int ind) {
    pos_2D = l.copy();
    index = ind;
    }
  void draw_vert(){
    noStroke();
    fill(255, 0, 0);
    translate(pos_2D.x,cylinder_height,pos_2D.y);
    rotateX(PI/2);
    circle(0,0,5);
    rotateX(-PI/2);
    translate(-pos_2D.x,-cylinder_height,-pos_2D.y);
    fill(255);
  }
}

class Edge implements Comparable<Edge>{
  Vertex start_pos;
  Vertex end_pos;
  float dist;
  float dist_to_goal;
  float cost;
  Edge(Vertex start,Vertex end) {
    start_pos = start;
    end_pos = end;
    dist = PVector.dist(start_pos.pos_2D, end_pos.pos_2D);
    dist_to_goal = PVector.dist(end_pos.pos_2D, endpoint.pos_2D);;
    }
  void draw_edge(){
    stroke(5);
    fill(255, 0, 0);
    translate(0,cylinder_height,0);
    rotateX(PI/2);
    line(start_pos.pos_2D.x, start_pos.pos_2D.y, end_pos.pos_2D.x, end_pos.pos_2D.y);
    rotateX(-PI/2);
    translate(0,-cylinder_height,0);
    fill(255);
    noStroke();
  }
  @Override     
  int compareTo(Edge e) {          
    return (this.cost < e.cost ? -1 : 
            (this.cost == e.cost ? 0 : 1));     
  }
}


//void add_edge(){
//  for (int i = 0; i < vertexs.size(); i++) {
//    Vertex point1 = vertexs.get(i);
//    TreeMap<Integer,Float> distance_map = new TreeMap<Integer,Float>(new Comparator<Integer>() {
//      @Override
//      public int compare(Integer o1, Integer o2) {
//        return o2.compareTo(o1);
//      }
//    });
//    for (int j = 0; j < vertexs.size(); i++) {
//      if (i!=j){
//        Vertex point2 = vertexs.get(j);
//        float dist = point2.pos_2D.dist(point1.pos_2D); 
//        distance_map.put(j,dist);
//      }
//    }
//    Iterator<Entry<Integer, Float>> it = distance_map.entrySet().iterator();
//    for (int k = 0; k < 5; i++) {
//      Entry<Integer, Float> entry = it.next();
//      int closest_key = entry.getKey();
//      Edge new_edge =new Edge(point1,vertexs.get(closest_key));
//      graph.addEdge(new_edge);
//    }
//  }
//}

int[] indexesOfTopElements(float[] orig, int nummax) {
  float[] copy = Arrays.copyOf(orig,orig.length);
  Arrays.sort(copy);
  float[] honey = Arrays.copyOfRange(copy,copy.length - nummax, copy.length);
  int[] result = new int[nummax];
  int resultPos = 0;
  for(int i = 0; i < orig.length; i++) {
    float onTrial = orig[i];
    int index = Arrays.binarySearch(honey,onTrial);
    if(index < 0) continue;
    result[resultPos++] = i;
  }
  return result;
}



ArrayList<Vertex> A_Star(Graph graph, int start_index, int end_index, int obj_num) {
  ArrayList<Vertex> solution = new ArrayList<Vertex>();
  ArrayList<Vertex> vertexs=graph.vertexs;
  ArrayList<Edge> edges=graph.edges;
  Queue<Vertex> queue = new LinkedList<Vertex>();
  //add start point to the queue
  queue.add(vertexs.get(start_index));
  vertexs.get(start_index).visited=true;
  float curr_cost = 0;
  while(!queue.isEmpty()) {
    Vertex vertex=queue.peek();
    ArrayList<Edge> valid_edge = new ArrayList<Edge>();
    for(Edge edge:edges){
      if(edge.start_pos.equals(vertex)&&edge.end_pos.visited==false){
        valid_edge.add(edge);
        edge.cost = curr_cost + edge.dist + edge.dist_to_goal;
      }
    }
    Collections.sort(valid_edge);
    for(Edge edge:valid_edge) {
      queue.add(edge.end_pos);
      edge.end_pos.visited=true;
      switch(obj_num){
        case 0:
          edge.end_pos.current_pre0 = edge.start_pos;
          if (edge.end_pos == vertexs.get(end_index)){
            Vertex vertex0 = vertex;
            print("find0!");
            while(vertex0.current_pre0 != null){
              solution.add(vertex0.current_pre0);
              vertex0 = vertex0.current_pre0;
            }
            return solution;
          }
        case 1:
          edge.end_pos.current_pre1 = edge.start_pos;
          if (edge.end_pos == vertexs.get(end_index)){
            Vertex vertex0 = vertex;
            print("find1!");
            while(vertex0.current_pre1 != null){
              solution.add(vertex0.current_pre1);
              vertex0 = vertex0.current_pre1;
            }
            return solution;
          }
          case 2:
          edge.end_pos.current_pre2 = edge.start_pos;
          if (edge.end_pos == vertexs.get(end_index)){
            Vertex vertex0 = vertex;
            print("find2!");
            while(vertex0.current_pre2 != null){
              solution.add(vertex0.current_pre2);
              vertex0 = vertex0.current_pre2;
            }
            return solution;
          }
      }
    }
    queue.poll();
  }
  return solution;
}

class Boid {
  PVector position;
  PVector velocity;
  PVector acceleration;
  int boid_color = color(255,255,255);
  float r;
  float maxforce;
  float maxspeed;
    Boid(float x, float y, PVector vel, int this_color) {
    acceleration = new PVector(0, 0);
    boid_color = this_color;
    float angle = random(PI);
    velocity = new PVector(0.05*vel.x+0.01*sin(angle), 0.05*vel.y+0.01*cos(angle));

    position = new PVector(x, y);
    r = 2.0;
    maxspeed = 1.0;
    maxforce = 0.03;
  }

  void run(ArrayList<Boid> boids) {
    flock(boids);
    update();
    borders();
    render();
  }
  
  void reset_velocity(PVector vel){
    velocity = vel;
  }

  void applyForce(PVector force) {
    acceleration.add(force);
  }

  void flock(ArrayList<Boid> boids) {
    PVector sep = separate(boids);
    PVector ali = align(boids);
    PVector coh = cohesion(boids);
    PVector col = avoid_collision();
    sep.mult(1.5);
    ali.mult(1.0);
    coh.mult(1.0);
    col.mult(1.0);
    applyForce(sep);
    applyForce(ali);
    applyForce(coh);
    applyForce(col);
  }
  void update() {
    velocity.add(acceleration);
    velocity.limit(maxspeed);
    position.add(velocity);
    acceleration.mult(0);
  }
  PVector seek(PVector target) {
    PVector desired = PVector.sub(target, position);
    // Scale to maximum speed
    desired.normalize();
    desired.mult(maxspeed);
    PVector steer = PVector.sub(desired, velocity);
    steer.limit(maxforce);
    return steer;
  }

  void render() {
    float theta = velocity.heading2D() + radians(90);
    //fill(200, 100);
    stroke(255);
    //noStroke();
    fill(boid_color);
    pushMatrix();
    translate(position.x, 20, position.y);
    rotate(theta);
    rotateX(PI/2);
    //rotateY(PI/2);
    beginShape(TRIANGLES);
    vertex(0, -r*2);
    vertex(-r, r*2);
    vertex(r, r*2);
    endShape();
    popMatrix();
  }


  void borders() {
    //if (position.x < -r) position.x = width+r;
    //if (position.y < -r) position.y = height+r;
    if (position.x > 200){
      velocity.x = 0;
      velocity.y = 0;
      acceleration.x = 0;
      acceleration.y = 0;
    }
  }
  PVector separate (ArrayList<Boid> boids) {
    float desiredseparation = 25.0f;
    PVector steer = new PVector(0, 0, 0);
    int count = 0;
    for (Boid other : boids) {
      float d = PVector.dist(position, other.position);
      if ((d > 0) && (d < desiredseparation)) {
        PVector diff = PVector.sub(position, other.position);
        diff.normalize();
        diff.div(d);
        steer.add(diff);
        count++;
      }
    }
    if (count > 0) {
      steer.div((float)count);
    }

    if (steer.mag() > 0) {
      steer.normalize();
      steer.mult(maxspeed);
      steer.sub(velocity);
      steer.limit(maxforce);
    }
    return steer;
  }
  PVector align (ArrayList<Boid> boids) {
    float neighbordist = 50;
    PVector sum = new PVector(0, 0);
    int count = 0;
    for (Boid other : boids) {
      float d = PVector.dist(position, other.position);
      if ((d > 0) && (d < neighbordist)) {
        sum.add(other.velocity);
        count++;
      }
    }
    if (count > 0) {
      sum.div((float)count);
      sum.normalize();
      sum.mult(maxspeed);
      PVector steer = PVector.sub(sum, velocity);
      steer.limit(maxforce);
      return steer;
    } 
    else {
      return new PVector(0, 0);
    }
  }
  PVector cohesion (ArrayList<Boid> boids) {
    float neighbordist = 50;
    PVector sum = new PVector(0, 0);
    int count = 0;
    for (Boid other : boids) {
      float d = PVector.dist(position, other.position);
      if ((d > 0) && (d < neighbordist)) {
        sum.add(other.position);
        count++;
      }
    }
    if (count > 0) {
      sum.div(count);
      return seek(sum);
    } 
    else {
      return new PVector(0, 0);
    }
  }
  
  PVector avoid_collision() {
    float dist_to_obs1 = PVector.dist(position, new PVector(0, 0))-20;
    float dist_to_obs2 = point_rec_dist(position,-100,0,50,100);
    PVector direct1 = PVector.sub(position, new PVector(0, 0));
    PVector direct2 = PVector.sub(position, new PVector(-100, 0));
    float f1 = 0;
    float f2 = 0;
    if (dist_to_obs1<5){
      f1= 1/dist_to_obs1;
    }
    if (dist_to_obs2<5){
      f2= 1/dist_to_obs1;
    }    
    return direct1.mult(f1).add(direct2.mult(f2));
  }
}


class Flock {
  ArrayList<Boid> boids;

  Flock() {
    boids = new ArrayList<Boid>();
  }
  
  void reset_boids_velocity(PVector vel){
    for (Boid b : boids) {
      b.reset_velocity(vel);
    }
  }

  void run() {
    for (Boid b : boids) {
      b.run(boids);
    }
  }

  void addBoid(Boid b) {
    boids.add(b);
  }

}
