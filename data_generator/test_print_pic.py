
# list = [[-0.016, 0.48, 0.032, 0.07], [-0.08, 0.48, 0.032, 0.07], [0.048, 0.48, 0.032, 0.07]]
#
# from tkinter import *
# from random import *
#
# tk = Tk()
# tk.title("apidemos.com")
# canvas = Canvas(tk, width=640, height=480)
# canvas.pack()
# for i in range(20):
#     x1, y1 = randint(1,640), randint(1,480)
#     x2, y2 = randint(1,640), randint(1,480)
#     if x1 > x2:
#         x1,x2 = x2,x1
#     if y1 > y2:
#         y1,y2 = y2,y1
#     canvas.create_rectangle(x1,y1,x2,y2)
#     # canvas.create_line(x1,y1,x2,y2)
#
# tk.mainloop()


# from tkinter import *
# from PIL import Image
# root = Tk()
# cv = Canvas(root)
# cv.create_rectangle(10,10,50,50)
# cv.pack()
# cv.postscript(file="a.eps",colormode='color')       #如将扩展名改为ps,出错信息相同
# img = Image.open("a.eps")
# img.save("a.png", "png")	#运行后，此句出错
# root.mainloop()
