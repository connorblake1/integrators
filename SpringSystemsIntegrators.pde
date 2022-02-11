ArrayList<PVector[]> x = new ArrayList<PVector[]>();
ArrayList<PVector[]> v = new ArrayList<PVector[]>();
PVector [][][] xm;
PVector [] vsm;
PFont font1;
PVector [] xsm;
PVector pivot;
int systems = 4;
float [] l = {200,200,200};
float [] k = {1,1,1};
float [] m = {1,1};
float uu = 0;
final float g = 200;
final float dt = .05;
final int cycle = 5;
int tick = 0;
final int s = 40;

void setup() {
  font1 = createFont("Arial", 30);
  textFont(font1);
  size(700,700);
  background(0);
  colorMode(HSB);
  pivot = new PVector(width/2,100);
  PVector sx0 = new PVector(width/2-100,400);
  PVector sx1 = new PVector(width/2+100,400);
  PVector sv0 = new PVector(-10,0);
  PVector sv1 = new PVector(50,35);
  uu  = getU(new PVector []{sx0,sx1}, new PVector[] {sv0,sv1});
  x.add(new PVector [] {new PVector(sx0.x,sx0.y),new PVector(sx1.x,sx1.y)});
  v.add(new PVector [] {new PVector(sv0.x,sv0.y),new PVector(sv1.x,sv1.y)});
  x.add(new PVector [] {new PVector(sx0.x,sx0.y),new PVector(sx1.x,sx1.y)});
  v.add(new PVector [] {new PVector(sv0.x,sv0.y),new PVector(sv1.x,sv1.y)});
  x.add(new PVector [] {new PVector(sx0.x,sx0.y),new PVector(sx1.x,sx1.y)});
  v.add(new PVector [] {new PVector(sx0.x-sv0.x*dt,sx0.y-sv0.y*dt),new PVector(sx1.x-sv1.x*dt,sx1.y-sv1.y*dt)});
  x.add(new PVector [] {new PVector(sx0.x,sx0.y),new PVector(sx1.x,sx1.y)});
  v.add(new PVector [] {new PVector(sv0.x,sv0.y),new PVector(sv1.x,sv1.y)});

}

void draw() {
  colorMode(HSB);
  fill(255,255,0);
  noStroke();
  rect(0,0,width,height);
  //fill(255);
  fill(64*0,255,255);
  text("Euler: " + round(100*(getU(x.get(0),v.get(0))-uu))/100.0,425,100);
  fill(64*1,255,255);
  text("Leapfrog: " + round(100*(getU(x.get(1),v.get(1))-uu))/100.0,425,150);
  fill(64*2,255,255);
  text("Verlet: " + round(100*(getUx(x.get(2),v.get(2))-uu))/100.0,425,50);
  fill(64*3,255,255);
  text("Yoshida: " + round(100.0*(getU(x.get(3),v.get(3))-uu))/100.0,425,200);

  Euler(v.get(0),x.get(0));
  LeapFrog(v.get(1),x.get(1));
  Verlet(v.get(2),x.get(2));
  Yoshida(v.get(3),x.get(3));
  drawSystemSkel(x.get(0),0);
  drawSystemSkel(x.get(1),1);
  drawSystemSkel(x.get(2),2);
  drawSystemSkel(x.get(3),3);
  
  //drawSystemFull(x.get(0));
}

PVector deepCopy(PVector in) {
  PVector r = new PVector(in.x,in.y,in.z);
  return r;}
  
float getU(PVector [] xin, PVector [] vin) {
    float u = 0;
    u+= sq(PVector.dist(xin[0],pivot)-l[0])*k[0]/2;
    u+= sq(PVector.dist(xin[1],pivot)-l[1])*k[1]/2;
    u+= sq(PVector.dist(xin[0],xin[1])-l[2])*k[2]/2;
    u += vin[0].magSq()/2*m[0];
    u += vin[1].magSq()/2*m[1];
    u -= m[0]*g*xin[0].y;
    u -= m[1]*g*xin[1].y;
  return u;
}
float getUx(PVector [] xin, PVector [] xino) {
    float u = 0;
    u+= sq(PVector.dist(xin[0],pivot)-l[0])*k[0]/2;
    u+= sq(PVector.dist(xin[1],pivot)-l[1])*k[1]/2;
    u+= sq(PVector.dist(xin[0],xin[1])-l[2])*k[2]/2;
    u += (PVector.mult(PVector.sub(xin[0],xino[0]),1/dt)).magSq()/2*m[0];
    u += (PVector.mult(PVector.sub(xin[0],xino[0]),1/dt)).magSq()/2*m[1];
    u -= m[0]*g*xin[0].y;
    u -= m[1]*g*xin[1].y;
  return u;
}

PVector[] calcF(PVector [] xin) {
  float dx0 = PVector.dist(xin[0],pivot)-l[0];
  float theta0 = atan2(pivot.y-xin[0].y, pivot.x-xin[0].x);
  float dx1 = PVector.dist(xin[1],pivot)-l[1];
  float theta1 = atan2(pivot.y-xin[1].y, pivot.x-xin[1].x);
  float dx2 = PVector.dist(xin[0],xin[1])-l[2];
  float theta2 = atan2(xin[0].y-xin[1].y,xin[0].x-xin[1].x);
  PVector [] f = new PVector[2];  
  f[0] = new PVector(cos(theta0)*k[0]*dx0-cos(theta2)*k[2]*dx2,sin(theta0)*k[0]*dx0-sin(theta2)*k[2]*dx2+m[0]*g);
  f[1] = new PVector(cos(theta1)*k[1]*dx1+cos(theta2)*k[2]*dx2,sin(theta1)*k[1]*dx1+sin(theta2)*k[2]*dx2+m[1]*g);
  return f;}
//system 0
void Euler(PVector [] vs, PVector [] xs) {
  for (int i = 0; i < 2; i ++) {
    vs[i] = PVector.add(vs[i], PVector.mult(calcF(xs)[i],dt));}
  for (int i = 0; i < 2; i++) {
    xs[i] = PVector.add(xs[i], PVector.mult(vs[i],dt));}
}
//system 1
void LeapFrog(PVector [] vs, PVector [] xs) {
  for (int i = 0; i < 2; i ++) {
    vs[i] = PVector.add(vs[i], PVector.mult(calcF(xs)[i],dt/2));}
  for (int i = 0; i < 2; i ++) {
    xs[i] = PVector.add(xs[i], PVector.mult(vs[i],dt));}
  for (int i = 0; i < 2; i ++) {
    vs[i] = PVector.add(vs[i], PVector.mult(calcF(xs)[i],dt/2));}
}
void Verlet(PVector [] xso, PVector [] xs) {
  PVector [] a = calcF(xs);
  deepCopySystem(2);
  for (int i = 0; i < 2; i ++) {
    xs[i] = new PVector(2*xsm[i].x-xso[i].x+dt*dt*a[i].x,2*xsm[i].y-xso[i].y+dt*dt*a[i].y);
  }
  xso[0] = new PVector(xsm[0].x,xsm[0].y);
  xso[1] = new PVector(xsm[1].x,xsm[1].y);
}

void Yoshida(PVector [] vs, PVector [] xs) {
  float w0 = -pow(2,.333)/(2-pow(2,.333));
  float w1 = 1/(2-pow(2,.333));
  float c1 = .6756;
  float c4 = c1;
  float c2 = -.1756;
  float c3 = c2;
  float d1 = w1;
  float d2 = w0;
  float d3 = d1;
  for (int i = 0; i < 2; i ++) {
    xs[i] = PVector.add(xs[i], PVector.mult(vs[i],c1*dt));}
  for (int i = 0; i < 2; i ++) {
    vs[i] = PVector.add(vs[i], PVector.mult(calcF(xs)[i],d1*dt));}
  for (int i = 0; i < 2; i ++) {
    xs[i] = PVector.add(xs[i], PVector.mult(vs[i],c2*dt));}
  for (int i = 0; i < 2; i ++) {
    vs[i] = PVector.add(vs[i], PVector.mult(calcF(xs)[i],d2*dt));}
  for (int i = 0; i < 2; i ++) {
    xs[i] = PVector.add(xs[i], PVector.mult(vs[i],c3*dt));}
  for (int i = 0; i < 2; i ++) {
    vs[i] = PVector.add(vs[i], PVector.mult(calcF(xs)[i],d3*dt));}
  for (int i = 0; i < 2; i ++) {
    xs[i] = PVector.add(xs[i], PVector.mult(vs[i],c4*dt));}
}

void deepCopySystem(int i) {
  xsm = new PVector[systems];
  vsm = new PVector[systems];
  for (int j = 0; j < 2; j++) {
    xsm[j] = deepCopy(x.get(i)[j]);
    vsm[j] = deepCopy(v.get(i)[j]);}
}

void drawSystemFull(PVector [] in) {
  colorMode(HSB);
  PVector [] ins = new PVector[]{in[0],pivot,in[1],pivot,in[0],in[1]};
  for (int i = 0; i < 3; i++) {
    PVector a = ins[i*2];
    PVector b = ins[i*2+1];
    float d = PVector.dist(a,b);
    resetMatrix();
    translate(a.x,a.y);
    rotate(atan2(b.y-a.y,b.x-a.x));
    for (float j = 0; j < 1; j+=.01) {
      stroke(100+100*sin(j*TWO_PI),255,255);
      line(map(j,0,1,0,d),-5,map(j,0,1,0,d),5);}}
  resetMatrix();
  colorMode(RGB);
  stroke(255);
  fill(255);
  ellipse(in[0].x,in[0].y,s,s);
  ellipse(in[1].x,in[1].y,s,s);
  fill(255,0,0);
  ellipse(pivot.x,pivot.y,s,s);}

void spring(PVector x1, PVector x2, float l) {
  translate(x1.x,x1.y);
  rotate(atan2(x2.y-x1.y,x2.x-x1.x));
  float d = PVector.dist(x1,x2);
  float waves = l/10;
  float amp = 5;
  for (int i = 0; i < d; i++) {
    point(i,amp*sin(i/d*TWO_PI*waves));}
  resetMatrix();
}
void drawSystemSkel(PVector [] in, int s) {
   //noStroke();
   //stroke(2*1.0/systems*255+tick*1.0/mem*255/systems,255,255);
   //line(pivot.x,pivot.y,in[0].x,in[0].y);
   //line(pivot.x,pivot.y,in[1].x,in[1].y);
   //line(in[0].x,in[0].y,in[1].x,in[1].y);
   colorMode(HSB);
   stroke(64*s,255,255);
   spring(pivot,in[0],l[0]);
   spring(pivot,in[1],l[1]);
   spring(in[0],in[1],l[2]);
   fill(64*s,255,255);
   noStroke();
   ellipse(in[0].x,in[0].y,20,20);
   ellipse(in[1].x,in[1].y,20,20);
   fill(255,0,255);
   ellipse(pivot.x,pivot.y,20,20);
}
