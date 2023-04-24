#contains function that does mce coord --> det coord conversion
#starting from Excel table: '150ghz_mapping.xlsx'
#SF - 2022/07


import numpy as np
import xlrd

xlrd.xlsx.ensure_elementtree_imported(False, None)
xlrd.xlsx.Element_has_iter = True

def mce2det(col, row): #works for col0-16

	excel_fn ='150ghz_mapping.xlsx'

	# Give the location of the file
	loc = (excel_fn)

	# To open Workbook
	wb = xlrd.open_workbook(loc)
	sheet = wb.sheet_by_index(0)

	# For row 0 and column 0
	#print(sheet.cell_value(1, 0)) #(row, col)

	#
	mcecol=[]
	mcerow=[]
	detcol=[]
	detrow=[]
	detpol=[]
	for i_row in range(2, sheet.nrows):
		mcecol.append(int(sheet.cell_value(i_row, 0)))
		mcerow.append(int(sheet.cell_value(i_row, 1)))
		try:
			detcol.append(int(sheet.cell_value(i_row, 2)))
			detrow.append(int(sheet.cell_value(i_row, 3)))
		except:
			detcol.append(sheet.cell_value(i_row, 2))
			detrow.append(sheet.cell_value(i_row, 3))
		detpol.append(sheet.cell_value(i_row, 4))

	nrows_det = len(detrow)
	ncol_det = len(detcol)

	mcerow=np.array(mcerow)
	mcecol=np.array(mcecol)

	detpol=np.array(detpol)
	detcol=np.array(detcol)
	detrow=np.array(detrow)

	print(mcerow)
	print(row)

	det_idx = np.where((mcerow==row) & (mcecol==col))

	drow = detrow[det_idx]
	dcol = detcol[det_idx]
	dpol = detpol[det_idx]

	print('number of dark det = ', len(detpol[np.where(detpol=='D')]))

	# tilemap = np.random.rand(nrows_det, ncol_det)
	# pl.imshow(tilemap)
	# pl.show()

	print('drow, dcol, dpol = ', drow, dcol, dpol)

	return dcol, drow, dpol
