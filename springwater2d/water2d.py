from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtGui import QVector2D
from OpenGL import GL
import math
class Mass:
	def __init__(self,mass,pos):
		self.mass=mass
		self.pos=pos
		self.vel=QVector2D(0,0)
		self.f=QVector2D(0,0)

		self.enter=False
	def Step(self,dt):
		#first we use forward euler
		a=self.f/self.mass
		self.vel=self.vel+a*dt
		self.vel.setX(0)
		self.pos=self.pos+self.vel*dt
		#reset for next step
		self.f=QVector2D(0,0)
	def OnEnterBox(self):
		if self.enter==True: return
		self.enter=True

	def OnExitBox(self):
		# print "exit"
		if self.enter==True:
			self.enter=False
			impulse = QVector2D(0,-100)
			self.ApplyImpulse(impulse)
	def ApplyImpulse(self,impulse):
		self.vel=self.vel+QVector2D(0,impulse.y())
class Spring:
	def __init__(self,k,damp,ma,mb):
		self.k=k
		self.damp=damp
		self.ma=ma
		self.mb=mb
		p0=self.ma.pos
		p1=self.mb.pos
		self.len=(p1-p0).length()
	def AddForce(self):
		#clac force and apply on mass point
		p0=self.ma.pos
		p1=self.mb.pos
		cur_len=(p1-p0).length()
		cur_dir=p1-p0
		cur_dir.normalize()
		dx=self.len-cur_len
		f=self.k*dx*cur_dir
		sign=1
		if self.ma.vel*cur_dir < 0: sign=-1
		self.ma.f=self.ma.f-f-self.ma.vel*cur_dir*self.damp*sign
		sign=1
		if self.mb.vel*cur_dir < 0: sign=-1
		self.mb.f=self.mb.f+f-self.mb.vel*cur_dir*self.damp*sign
	def DebugDraw(self,scale):
		p0=self.ma.pos*scale
		p1=self.mb.pos*scale
		GL.glBegin(GL.GL_LINES)
		GL.glVertex2d(p0.x(),p0.y())
		GL.glVertex2d(p1.x(),p1.y())
		GL.glEnd()
class Box:
	width=10
	height=10
	def __init__(self,pos,f):
		self.pos=pos
		self.f=f
		self.vel=QVector2D(0,0)
	def Step(self,dt):
		# euler integrate
		self.vel=self.vel+self.f*dt
		self.pos=self.pos+self.vel*dt
	def Left(self):
		return self.pos.x()-self.width/2
	def Right(self):
		return self.pos.x()+self.width/2
	def Bottom(self):
		return self.pos.y()-self.height/2
	def Top(self):
		return self.pos.y()+self.height/2
	def PointInside(self,p):
		if p.x()>self.Left() and p.x()<self.Right():
			if p.y()>self.Bottom() and p.y()<self.Top():
				return True
		return False
	def Draw(self):
		l=self.Left()
		r=self.Right()
		b=self.Bottom()
		t=self.Top()
		GL.glVertex2f(l,b)
		GL.glVertex2f(r,b)
		GL.glVertex2f(r,t)
		GL.glVertex2f(l,t)

class Water2d:
	mass_val = 0.01
	k		 = 10
	damp	 = 0.05
	iter_num = 8
	spread	 = 0.3
	gravity	 = -980
	impause_factor = 1
	# box just for test
	def __init__(self,width,height,grid_num):
		self.width=width
		self.height=height
		self.ratio=height/float(width)
		self.grid_num=grid_num
		self.grid_step=100.0/grid_num
		self.masses=[]
		self.static_masses=[]
		self.springs=[]
		# boxes 
		self.boxes=[]
		for i in range(grid_num+1):
			mass = Mass(self.mass_val,QVector2D(self.grid_step*i,100*self.ratio))
			self.masses.append(mass)
			mass = Mass(self.mass_val,QVector2D(self.grid_step*i,0))
			self.static_masses.append(mass)
		for i in range(grid_num):
			# self.springs.append(Spring(self.k,0,self.masses[i]       ,self.masses[i+1]))
			self.springs.append(Spring(self.k,self.damp,self.static_masses[i],self.masses[i]))
			# self.springs.append(Spring(self.k,self.static_masses[i],self.masses[i+1]))
	def Update(self,dt):
		#calc spring scale
		for spring in self.springs:
			spring.AddForce()
		self.__collision_detect()
		for i in range(0,self.grid_num+1):
			self.masses[i].Step(dt)
		# convolution
		dl=[0 for x in range(self.grid_num+1)] 
		dr=[0 for x in range(self.grid_num+1)]
		for it in range(self.iter_num):
			for i in range (self.grid_num+1):
				if i>0:
					dl[i]=self.spread*(self.masses[i].pos.y()-self.masses[i-1].pos.y())
					self.masses[i-1].vel=self.masses[i-1].vel+QVector2D(0,dl[i])
				if i<self.grid_num:
					dr[i]=self.spread*(self.masses[i].pos.y()-self.masses[i+1].pos.y())
					self.masses[i+1].vel=self.masses[i+1].vel+QVector2D(0,dr[i])
			for i in range (self.grid_num):
				if i>0:
					self.masses[i-1].pos=self.masses[i-1].pos+QVector2D(0,dl[i])
				if i<self.grid_num:
					self.masses[i+1].pos=self.masses[i+1].pos+QVector2D(0,dr[i])
		# --------------------------------------
		for b in self.boxes:
			b.Step(dt)
	def AddBox(self,pos):
		box=Box(pos,QVector2D(0,self.gravity))
		box.vel=QVector2D(0,-50)
		self.boxes.append(box)
	def __collision_detect(self):
		for b in self.boxes:
			l=int(b.Left()/self.grid_step)+1
			r=int(b.Right()/self.grid_step)-1
			# print "l: "+str(l)
			# print "r: "+str(r)
			for i in range(l,r):
				m=self.masses[i]
				if b.PointInside(m.pos):
					m.OnEnterBox()
				else:
					m.OnExitBox()
	def Render(self):
		scale = self.width/100.0
		GL.glScalef(scale,scale,1)
		GL.glBegin(GL.GL_LINES)
		
		for i in range(self.grid_num):
			p0=self.masses[i].pos
			p1=self.masses[i+1].pos
			GL.glVertex2f(p0.x(),p0.y())
			GL.glVertex2f(p1.x(),p1.y())
		GL.glEnd()
		GL.glBegin(GL.GL_QUADS)
		for b in self.boxes:
			b.Draw()
		GL.glEnd()
		GL.glScalef(1,1,1)
	def HandleMouseMove(self,x,y,dx,dy):
		simx=x*100/self.width
		mid=int(simx%self.grid_step)
		mid = self.grid_num/2
		if mid<=0: mid=1
		if mid>=self.grid_num: mid=self.grid_num-1
		for i in range(max(mid-5,1),min(mid+5,self.grid_num-1)):
			self.masses[i].ApplyImpulse(QVector2D(dx,dy))