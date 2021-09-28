import matplotlib.pyplot as plt
import matplotlib
from matplotlib import gridspec
from matplotlib import rc
#rc('text', usetex=True)
import jjrplot as jplt
import numpy as np

def img(x1d, x2d, Vals, bSaveFig, bPNG, filename, cmapstr,
		xlabel, ylabel, title_str, cblabelstr, vminV, vmaxV, t):

	#cmapstr = 'viridis'

	#fig, ax = plt.subplots(figsize=(24,6))
	fig, ax = plt.subplots()

	imgplot = plt.imshow((Vals),cmap=cmapstr,
						extent = (x1d.min(), x1d.max(),
						x2d.min(), x2d.max()),
						aspect='auto',
						origin='lower',
						vmin=vminV, vmax=vmaxV)

	ax = jplt.axAdj(ax)


	plt.xlabel(r'$' + xlabel + '$', fontsize=20)
	plt.ylabel(r'$' + ylabel + '$', fontsize=20)
	plt.title(title_str, fontsize=16, loc='left')

	trans = ax.get_yaxis_transform() # y in data untis, x in axes fraction
	plt.annotate(r'$t(\Omega^{{-1}})={:.2f}$'.format(t),xy=(0.65,x2d.max()*1.05), 
						xycoords=trans, fontsize=20)
	
	cb_label = r'' + cblabelstr

	bSetTicks = True
	if(bSetTicks):
		cb = plt.colorbar()
	else:
		cb = plt.colorbar()
		ax.set_yticklabels([])
		ax.set_xticklabels([])

	fg_color = "black"
	# set colorbar label plus label color
	cb.set_label(cb_label, color=fg_color, fontsize=20)

	# set colorbar tick color
	cb.ax.yaxis.set_tick_params(color=fg_color, labelsize=16)

	# set colorbar edgecolor 
	cb.outline.set_edgecolor(fg_color)

	# set colorbar ticklabels
	plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), color=fg_color)

	if(bSaveFig == 0):
		plt.show()
	if(bSaveFig == 1):
		if(bPNG):
			picstr = filename + '.png'
			plt.savefig(picstr, bbox_inches='tight')
		else:
			picstr = filename + '.pdf'
			plt.savefig(picstr, bbox_inches='tight')
			# can add dpi=600 , 800, here

	plt.clf()
	plt.close()


def twoimg(x1d, x2d, Vals1, Vals2, bSaveFig, bPNG, filename, cmapstr,
		xlabel, ylabel, title_str, cblabelstr, vminV, vmaxV, t):

	fig, [ax1, ax2] = plt.subplots(1,2)#,sharey=True)

	imgplot1 = ax1.imshow((Vals1),cmap=cmapstr,
		extent = (x1d.min(), x1d.max(),
		x2d.min(), x2d.max()),
		#extent = (-0.1, 0.1,
		#-0.1, 0.1),
		origin='lower',
		vmin=vminV, vmax=vmaxV)

	ax1 = jplt.axAdj(ax1)

	ax1.set_xlabel(r'$' + xlabel + '$', fontsize=20)
	ax1.set_ylabel(r'$' + ylabel + '$', fontsize=20)
	plt.title(title_str, fontsize=16, loc='left')


	imgplot2 = ax2.imshow((Vals2),cmap=cmapstr,
		extent = (x1d.min(), x1d.max(),
		x2d.min(), x2d.max()),
		origin='lower',
		vmin=vminV, vmax=vmaxV)

	ax2 = jplt.axAdj(ax2)

	ax2.set_xlabel(r'$' + xlabel + '$', fontsize=20)
	ax2.set_yticklabels([])

	trans = ax1.get_yaxis_transform() # y in data untis, x in axes fraction
	plt.annotate(r'\textrm{no self grav}',xy=(0.01,x2d.max()*1.05), 
						xycoords=trans, fontsize=20)

	trans = ax2.get_yaxis_transform() # y in data untis, x in axes fraction
	plt.annotate(r'$t(\Omega^{{-1}})={:.2f}$'.format(t),xy=(-0.3,x2d.max()*1.05), 
						xycoords=trans, fontsize=20)
	plt.annotate(r'\textrm{self grav}',xy=(0.60,x2d.max()*1.05), 
						xycoords=trans, fontsize=20)

	cb_label = r'' + cblabelstr

	#cbaxes = fig.add_axes([0.8, 0.1, 0.03, 0.8])

	fig.subplots_adjust(wspace=0.05)

	fig.subplots_adjust(top=0.92)
	cbar_ax = fig.add_axes([0.12, 0.84, 0.79, 0.02])

	#fig.subplots_adjust(right=0.8)
	#cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])

	bSetTicks = True
	if(bSetTicks):
		cb = fig.colorbar(imgplot2, orientation='horizontal',
			cax = cbar_ax)
	else:
		cb = plt.colorbar()
		ax2.set_yticklabels([])
		ax2.set_xticklabels([])
		ax1.set_yticklabels([])
		ax1.set_xticklabels([])

	fg_color = "black"
	# set colorbar label plus label color
	cb.set_label(cb_label, color=fg_color, fontsize=20)

	cb.ax.xaxis.set_ticks_position('top')
	cb.ax.xaxis.set_label_position('top')

	# set colorbar tick color
	cb.ax.xaxis.set_tick_params(color=fg_color, labelsize=16)

	# set colorbar edgecolor 
	cb.outline.set_edgecolor(fg_color)

	# set colorbar ticklabels
	plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), color=fg_color)

	if(bSaveFig == 0):
		plt.show()
	if(bSaveFig == 1):
		if(bPNG):
			picstr = filename + '.png'
			plt.savefig(picstr, bbox_inches='tight')
		else:
			picstr = filename + '.pdf'
			plt.savefig(picstr, dpi=800, bbox_inches='tight')

	plt.clf()
	plt.close()

def twoimgvert(x1d, x2d1, x2d2, Vals1, Vals2, bSaveFig, bPNG, filename, cmapstr,
		xlabel, ylabel, title_str, cblabelstr, vminV, vmaxV, t):

	fig = plt.figure()#(2,1)#sharex=True)
	gs = gridspec.GridSpec(2, 1, height_ratios=[0.25, 1], figure=fig)
	ax1 = plt.subplot(gs[0])
	ax2 = plt.subplot(gs[1])

	imgplot1 = ax1.imshow((Vals1),cmap=cmapstr,
		extent = (x1d.min(), x1d.max(),
		x2d1.min(), x2d1.max()),
		aspect='auto',
		#extent5= (-0.1, 0.1,
		#-0.1, 0.1),
		origin='lower',
		vmin=vminV, vmax=vmaxV)

	#ax1 = plt.gca()
	yticks = [-0.1,0,0.1]
	ax1.set_yticks(yticks)

	ax1 = jplt.axAdj(ax1)
	ax1.set_xticklabels([])

	ax1.set_ylabel(r'$' + ylabel + '$', fontsize=20)
	plt.title(title_str, fontsize=16, loc='left')


	imgplot2 = ax2.imshow((Vals2),cmap=cmapstr,
		extent = (x1d.min(), x1d.max(),
		x2d2.min(), x2d2.max()),
		aspect='auto',
		origin='lower',
		vmin=vminV, vmax=vmaxV)

	ax2 = jplt.axAdj(ax2)

	ax2.set_ylabel(r'$' + ylabel + '$', fontsize=20)
	ax2.set_xlabel(r'$' + xlabel + '$', fontsize=20)

	cb_label = r'' + cblabelstr

	#cbaxes = fig.add_axes([0.8, 0.1, 0.03, 0.8])

	fig.subplots_adjust(wspace=0.0)

	fig.subplots_adjust(right=0.8)
	cbar_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])

	#fig.subplots_adjust(right=0.8)
	#cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])

	bSetTicks = True
	if(bSetTicks):
		cb = fig.colorbar(imgplot2,cax = cbar_ax)
	else:
		cb = plt.colorbar()
		ax2.set_yticklabels([])
		ax2.set_xticklabels([])
		ax1.set_yticklabels([])
		ax1.set_xticklabels([])

	fg_color = "black"
	# set colorbar label plus label color
	cb.set_label(cb_label, color=fg_color, fontsize=20)

	#cb.ax.xaxis.set_ticks_position('top')
	#cb.ax.xaxis.set_label_position('top')

	# set colorbar tick color
	cb.ax.yaxis.set_tick_params(color=fg_color, labelsize=18)

	# set colorbar edgecolor 
	cb.outline.set_edgecolor(fg_color)

	# set colorbar ticklabels
	plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), color=fg_color)

	if(bSaveFig == 0):
		plt.show()
	if(bSaveFig == 1):
		if(bPNG):
			picstr = filename + '.png'
			plt.savefig(picstr, bbox_inches='tight')
		else:
			picstr = filename + '.pdf'
			plt.savefig(picstr, bbox_inches='tight')

	plt.clf()
	plt.close()

def threeimgvert(x1d, x2d1, x2d2, x2d3, Vals1, Vals2, Vals3, 
		bSaveFig, bPNG, filename, cmapstr,
		xlabel, ylabel, title_str, cblabelstr, vminV, vmaxV, t):

	fig = plt.figure(figsize=[6,9])#(2,1)#sharex=True)
	gs = gridspec.GridSpec(3, 1, height_ratios=[0.25, 0.5, 1], figure=fig)
	ax1 = plt.subplot(gs[0])
	ax2 = plt.subplot(gs[1])
	ax3 = plt.subplot(gs[2])


# ------------------ 1st Plot --------------------------------------------- #

	imgplot1 = ax1.imshow((Vals1),cmap=cmapstr,
		extent = (x1d.min(), x1d.max(),
		x2d1.min(), x2d1.max()),
		aspect='auto',
		#extent5= (-0.1, 0.1,
		#-0.1, 0.1),
		origin='lower',
		vmin=vminV, vmax=vmaxV)

	#ax1 = plt.gca()
	yticks = [-0.1,0.1]
	ax1.set_yticks(yticks)

	ax1 = jplt.axAdj(ax1)
	ax1.set_xticklabels([])

	ax1.set_ylabel(r'$' + ylabel + '$', fontsize=20)


# ------------------ 2nd Plot --------------------------------------------- #

	imgplot2 = ax2.imshow((Vals2),cmap=cmapstr,
		extent = (x1d.min(), x1d.max(),
		x2d2.min(), x2d2.max()),
		aspect='auto',
		origin='lower',
		vmin=vminV, vmax=vmaxV)
	yticks = [-0.2,0,0.2]
	ax2.set_yticks(yticks)

	ax2 = jplt.axAdj(ax2)
	ax2.set_xticklabels([])

	ax2.set_ylabel(r'$' + ylabel + '$', fontsize=20)

# ------------------ 3rd Plot --------------------------------------------- #


	imgplot3 = ax3.imshow((Vals3),cmap=cmapstr,
		extent = (x1d.min(), x1d.max(),
		x2d3.min(), x2d3.max()),
		aspect='auto',
		origin='lower',
		vmin=vminV, vmax=vmaxV)
	yticks = [-0.4,-0.2,0,0.2,0.4]
	ax3.set_yticks(yticks)

	ax3 = jplt.axAdj(ax3)

	ax3.set_ylabel(r'$' + ylabel + '$', fontsize=20)
	ax3.set_xlabel(r'$' + xlabel + '$', fontsize=20)

# ----------------------------------------------------------------------- #

	plt.title(title_str, fontsize=16, loc='left')

	cb_label = r'' + cblabelstr

	#cbaxes = fig.add_axes([0.8, 0.1, 0.03, 0.8])

	fig.subplots_adjust(wspace=0.0)

	fig.subplots_adjust(right=0.8)
	cbar_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])

	#fig.subplots_adjust(right=0.8)
	#cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])

	bSetTicks = True
	if(bSetTicks):
		cb = fig.colorbar(imgplot2,cax = cbar_ax)
	else:
		cb = plt.colorbar()
		ax2.set_yticklabels([])
		ax2.set_xticklabels([])
		ax1.set_yticklabels([])
		ax1.set_xticklabels([])

	fg_color = "black"
	# set colorbar label plus label color
	cb.set_label(cb_label, color=fg_color, fontsize=20)

	#cb.ax.xaxis.set_ticks_position('top')
	#cb.ax.xaxis.set_label_position('top')

	# set colorbar tick color
	cb.ax.yaxis.set_tick_params(color=fg_color, labelsize=18)

	# set colorbar edgecolor 
	cb.outline.set_edgecolor(fg_color)

	# set colorbar ticklabels
	plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), color=fg_color)

	if(bSaveFig == 0):
		plt.show()
	if(bSaveFig == 1):
		if(bPNG):
			picstr = filename + '.png'
			plt.savefig(picstr, dpi=1200, bbox_inches='tight')
		else:
			picstr = filename + '.pdf'
			plt.savefig(picstr, dpi=1200, bbox_inches='tight')

	plt.clf()
	plt.close()


def img_wcircs(x1d, x2d, Vals, 
		xcirc, ycirc,
		bSaveFig, bPNG, filename, cmapstr,
		xlabel, ylabel, title_str, cblabelstr, vminV, vmaxV, t):

	#fig, ax = plt.subplots(figsize=(6,6))
	fig, ax = plt.subplots()

	ax.scatter(xcirc, ycirc, s=190*(0.00388205/0.00388205)**2, lw=2, facecolors='none', edgecolors='w')

	imgplot = ax.imshow((Vals),cmap=cmapstr,
						extent = (x1d.min(), x1d.max(),
						x2d.min(), x2d.max()),
						aspect='auto',
						origin='lower',
						vmin=vminV, vmax=vmaxV)


	ax = jplt.axAdj(ax)


	plt.xlabel(r'$' + xlabel + '$', fontsize=20)
	plt.ylabel(r'$' + ylabel + '$', fontsize=20)
	plt.title(title_str, fontsize=16, loc='left')

	trans = ax.get_yaxis_transform() # y in data untis, x in axes fraction
	plt.annotate(r'$t(\Omega^{{-1}})={:.2f}$'.format(t),xy=(0.65,x2d.max()*1.05), 
						xycoords=trans, fontsize=20)
	
	cb_label = r'' + cblabelstr

	bSetTicks = True
	if(bSetTicks):
		cb = plt.colorbar(imgplot)
	else:
		cb = plt.colorbar()
		ax.set_yticklabels([])
		ax.set_xticklabels([])

	fg_color = "black"
	# set colorbar label plus label color
	cb.set_label(cb_label, color=fg_color, fontsize=20)

	# set colorbar tick color
	cb.ax.yaxis.set_tick_params(color=fg_color, labelsize=16)

	# set colorbar edgecolor 
	cb.outline.set_edgecolor(fg_color)

	# set colorbar ticklabels
	plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), color=fg_color)

	if(bSaveFig == 0):
		plt.show()
	if(bSaveFig == 1):
		if(bPNG):
			picstr = filename + '.png'
			plt.savefig(picstr, dpi=800, bbox_inches='tight')
		else:
			picstr = filename + '.pdf'
			plt.savefig(picstr, dpi=800, bbox_inches='tight')

	plt.clf()
	plt.close()

def img_wcircs_warrow(x1d, x2d, Vals, 
		xcirc, ycirc, RH, vxcirc, vycirc,
		bSaveFig, bPNG, filename, cmapstr,
		xlabel, ylabel, title_str, cblabelstr, vminV, vmaxV, t):

	#fig, ax = plt.subplots(figsize=(7,6))
	fig, ax = plt.subplots()

	#ax.scatter(xcirc, ycirc, s=190*(RH/0.00388205)**2, lw=2, facecolors='none', edgecolors='w')

	if(len(vxcirc) > 0):
		Larrow = (x2d.max() - x2d.min())/12.
		vplt = np.sqrt(vxcirc[0]**2 + vycirc[0]**2)

		for i in range(len(vxcirc)):

			dx = Larrow*(vxcirc[i]/vplt)
			dy = Larrow*(vycirc[i]/vplt)

			#ax.arrow(xcirc[i], ycirc[i], dx, dy, ec=None, fc='k')
	#ax.plot([xcirc, xcirc+RH], [ycirc, ycirc], 'k-', lw=2)

	imgplot = ax.imshow((Vals),cmap=cmapstr,
						extent = (x1d.min(), x1d.max(),
						x2d.min(), x2d.max()),
						aspect='auto',
						origin='lower',
						vmin=vminV, vmax=vmaxV)

	circles = [plt.Circle((xi,yi), radius=ri) for xi,yi,ri in zip(xcirc,ycirc,RH)]
	c = matplotlib.collections.PatchCollection(circles, edgecolor='w', facecolor='', lw=1)
	ax.add_collection(c)

	#fig.subplots_adjust(left=0,right=1,bottom=0,top=1)

	ax.set_aspect('equal','box')

	ax = jplt.axAdj(ax)

	plt.xlabel(r'$' + xlabel + '$', fontsize=20)
	plt.ylabel(r'$' + ylabel + '$', fontsize=20)
	plt.title(title_str, fontsize=16, loc='left')

	trans = ax.get_yaxis_transform() # y in data untis, x in axes fraction
	plt.annotate(r'$t(\Omega^{{-1}})={:.2f}$'.format(t),xy=(0.65,x2d.max()*1.05), 
						xycoords=trans, fontsize=20)
	
	cb_label = r'' + cblabelstr

	bSetTicks = True
	if(bSetTicks):
		cb = plt.colorbar(imgplot)
	else:
		cb = plt.colorbar()
		ax.set_yticklabels([])
		ax.set_xticklabels([])

	fg_color = "black"
	# set colorbar label plus label color
	cb.set_label(cb_label, color=fg_color, fontsize=20)

	# set colorbar tick color
	cb.ax.yaxis.set_tick_params(color=fg_color, labelsize=16)

	# set colorbar edgecolor 
	cb.outline.set_edgecolor(fg_color)

	# set colorbar ticklabels
	plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), color=fg_color)

	#ax.axis("equal")

	if(bSaveFig == 0):
		plt.show()
	if(bSaveFig == 1):
		if(bPNG):
			picstr = filename + '.png'
			plt.savefig(picstr, dpi=800, bbox_inches='tight')
		else:
			picstr = filename + '.pdf'
			plt.savefig(picstr, dpi=800, bbox_inches='tight')

	plt.clf()
	plt.close()

def img_wpars(x1d, x2d, Vals, px1, px2, px3, bSaveFig, filename,
						xlabel, ylabel, cblabelstr, vminV, vmaxV, t):

	fig, ax = plt.subplots(figsize=(6,6))

	imgplot = plt.imshow((Vals),extent = (x1d.min(), x1d.max(),
						x2d.min(), x2d.max()),
						origin='lower',
						vmin=vminV, vmax=vmaxV)

	ax.plot(px1[np.where(px3>0.0)], px2[np.where(px3>0.0)], 'k.', ms=1)

	ax = jplt.axAdj(ax)


	plt.xlabel(r'$' + xlabel + '$', fontsize=20)
	plt.ylabel(r'$' + ylabel + '$', fontsize=20)

	plt.annotate(r'$t={:.2f}$'.format(t),(0.0,x2d.max()*0.999), fontsize=20)
	
	cb_label = r'' + cblabelstr

	cb = plt.colorbar()
	fg_color = "black"
	# set colorbar label plus label color
	cb.set_label(cb_label, color=fg_color, fontsize=20)

	# set colorbar tick color
	cb.ax.yaxis.set_tick_params(color=fg_color, labelsize=16)

	# set colorbar edgecolor 
	cb.outline.set_edgecolor(fg_color)

	# set colorbar ticklabels
	plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), color=fg_color)

	if(bSaveFig == 0):
		plt.show()
	if(bSaveFig == 1):
		picstr = filename + '.png'
		plt.savefig(picstr, bbox_inches='tight')

	plt.clf()
	plt.close()
