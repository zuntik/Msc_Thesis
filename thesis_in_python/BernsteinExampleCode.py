from bernsteinlib import *
import matplotlib.pyplot as plt

T = 3
t = np.linspace(0, T, 100)
p = np.array([2, 4, -1, 2, 5]).reshape((-1, 1))

_, ax = plt.subplots()
plt.title('First example')
# ax.xticks(rotation=90)
ax.invert_yaxis()
ax.set_xlabel('t')
ax.set_ylabel('b(t)')
ax.axis('equal')
plt.figtext(0.8, 0.8, "the first example")

curveplot, ctrlpoints = bernsteinPlot(p, T, ax=ax)

curveplot.set_label('a')
ctrlpoints.set_label('b')
ax.legend()

plt.show()


p2 = np.array([[0, np.cos(1), np.cos(2), np.cos(3), np.cos(4), np.cos(5), 0],
               [0, np.sin(1), np.sin(2), np.sin(3), np.sin(4), np.sin(5), 0]]).T

bernsteinPlot(p2, 2)


fig, ax = plt.subplots()

pointsa = bernsteinEval(p, T, t)

ax.plot(t, pointsa ** 2 - 1)
bernsteinPlot(bernsteinPow(p, 2), T, ax=ax)

# Test to mon
_, ax = plt.subplots()
m = bernsteinToMon(p, T)

ax.plot(t, bernsteinEval(p - 0.1, T, t))
ax.plot(t, [np.polyval(m, i) for i in t])


# Test Derivative

plt.figure()

plt.plot(t, bernsteinEval(bernsteinDeriv(p, T) - 1, T, t))
plt.plot(t, [np.polyval(np.polyder(m.flatten()), i) for i in t])


# Test Integral

plt.figure()
m_int = np.polyint(m.flatten())
p_int = bernsteinAntiDeriv(p, T, 0)

plt.plot(t, bernsteinEval(p_int - 0.2, T, t))
plt.plot(t, [np.polyval(m_int, i) for i in t])


# Test 3-D

p3 = np.array([[0, 0, 0], [1, np.cos(1), np.sin(1)], [2, np.cos(2), np.sin(2)], [3, np.cos(3), np.sin(3)],
               [4, np.cos(4), np.sin(4)], [5, np.cos(5), np.sin(5)], [6, np.cos(6), np.sin(6)], [7, 0, 0]])

bernsteinPlot(p3, 1)

plt.show()

bernsteinDegrElevMat(5, 6)
