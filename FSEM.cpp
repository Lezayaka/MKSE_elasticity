#include "Headers.h"

Point zero(const Point& p) {
	return { 0, 0 };
}


//FSEM::FSEM(double E, double nu, const Point& a, const Point& b, 
//	size_t n_x, size_t n_y, int coef_val_x, int coef_val_y): 
//	E(E), nu(nu), a(a), b(b), n_side_x(n_x - 1), 
//	n_side_y(n_y - 1), coef_x(coef_val_x), coef_y(coef_val_y),
//	fem(a, b, (n_y - 1)* coef_y + 1, (n_x - 1)* coef_x + 1) {
FSEM::FSEM(double E, double nu, const Point& a, const Point& b,
	size_t n_x, size_t n_y, int coef_x, int coef_y) :
	a(a), b(b), E(E), nu(nu), n_side_x(n_x - 1),
	n_side_y(n_y - 1), 
	fem(a, b, (n_y - 1)* coef_y + 1, (n_x - 1)* coef_x + 1),
	K((coef_x * (n_x - 1) + 1) * (coef_y * (n_y - 1) + 1)),
	f((coef_x * (n_x - 1) + 1) * (coef_y* (n_y - 1) + 1))
	{

	double h_x = (b.x - a.x) / n_side_x, h_y = (b.y - a.y) / n_side_y; // шаг
	
	int n_nodes = 2 * (n_side_x + n_side_y);
	nodes.resize(n_nodes);
	
	basis = std::vector<std::vector<Point>>(2 * n_nodes,
		std::vector<Point>((coef_x * n_side_x + 1) * (coef_y * n_side_y + 1)));

	basis_coefficients.resize(2 * n_nodes);
	
	for (int i = 0; i < n_side_y; ++i) {
		nodes[i] = { a.x, a.y + i * h_y };
		nodes[n_side_x + n_side_y + i] = { b.x, b.y - i * h_y };
	}
	for (int i = 0; i < n_side_x; ++i) {
		nodes[n_side_y + i] = { a.x + i * h_x, b.y };
		nodes[n_side_x + 2 * n_side_y + i] = { b.x - i * h_x, a.y };
	}

}

Point& FSEM::get_node(size_t n) {
	return nodes[n];
}

const Point& FSEM::get_node(size_t n) const {
	return nodes[n];
}

void FSEM::print_nodes() const {
	for (size_t i = 0; i < nodes.size(); ++i) {
		std::cout << (*this)[i] << "\t";
	}
	std::cout << '\n';
}

void FSEM::construct_basis() {
	// количество узлов сетки для МКЭ
	/*int n_x = n_side_x * coef_x + 1, n_y = n_side_y * coef_y + 1;
	double dx = (b.x - a.x) / (n_x - 1), dy = (b.y - a.y) / (n_y - 1);*/
	double dx = (b.x - a.x) / n_side_x, dy = (b.y - a.y) / n_side_y; // шаг

	//******************ЛЕВЫЙ НИЖНИЙ УГОЛ******************
	
	// находим базисную функцию для, равную 1 в a по х компоненте
	fem.set_boundaries('W',
		[&](const Point& p)
		{
			double neighbor_y = a.y + dy;
			if (p.y < neighbor_y)
				return Point{ (neighbor_y - p.y) / dy, 0 };

			return Point{ 0, 0 };
		});
	fem.set_boundaries('N', zero);
	fem.set_boundaries('E', zero);
	fem.set_boundaries('S',
		[&](const Point& p)
		{
			double neighbor_x = a.x + dx;
			if (p.x < neighbor_x)
				return Point{ (neighbor_x - p.x) / dx, 0 };

			return Point{ 0, 0 };
		});

	fem.construct_AF(E, nu, zero);
	auto [A, F] = fem.get_AF();
	fem.apply_boundaries();
	
	basis[0] = fem.solve();
	
	// находим базисную функцию для, равную 1 в a по y компоненте
	fem.set_boundaries('W',
		[&](const Point& p)
		{
			double neighbor_y = a.y + dy;
			if (p.y < neighbor_y)
				return Point{ 0, (neighbor_y - p.y) / dy};

			return Point{ 0, 0 };
		});
	fem.set_boundaries('N', zero);
	fem.set_boundaries('E', zero);
	fem.set_boundaries('S',
		[&](const Point& p)
		{
			double neighbor_x = a.x + dx;
			if (p.x < neighbor_x)
				return Point{ 0, (neighbor_x - p.x) / dx };

			return Point{ 0, 0 };
		});
	fem.set_AF(A, F);
	fem.apply_boundaries();

	basis[1] = fem.solve();
	
	//******************ЛЕВЫЙ ВЕРХНИЙ УГОЛ******************
	
	// находим базисную функцию для, равную 1 в {a.x, b.y} по х компоненте
	fem.set_boundaries('W',
		[&](const Point& p)
		{
			double neighbor_y = b.y - dy;
			if (p.y > neighbor_y)
				return Point{ (p.y - neighbor_y) / dy, 0 };

			return Point{ 0, 0 };
		});
	fem.set_boundaries('N',
		[&](const Point& p)
		{
			double neighbor_x = a.x + dx; 
			if (p.x < neighbor_x)
				return Point{ (neighbor_x - p.x) / dx, 0 };

			return Point{ 0, 0 };
		});
	fem.set_boundaries('E', zero);
	fem.set_boundaries('S', zero);

	fem.set_AF(A, F);
	fem.apply_boundaries();
	basis[2 * n_side_y] = fem.solve();

	// находим базисную функцию для, равную 1 в {a.x, b.y} по y компоненте
	fem.set_boundaries('W',
		[&](const Point& p)
		{
			double neighbor_y = b.y - dy;
			if (p.y > neighbor_y)
				return Point{ 0, (p.y - neighbor_y) / dy};

			return Point{ 0, 0 };
		});
	fem.set_boundaries('N',
		[&](const Point& p)
		{
			double neighbor_x = a.x + dx;
			if (p.x < neighbor_x)
				return Point{ 0, (neighbor_x - p.x) / dx};

			return Point{ 0, 0 };
		});
	fem.set_boundaries('E', zero);
	fem.set_boundaries('S', zero);

	fem.set_AF(A, F);
	fem.apply_boundaries();
	basis[2 * n_side_y + 1] = fem.solve();

	//******************ПРАВЫЙ ВЕРХНИЙ УГОЛ******************
	
	// находим базисную функцию для, равную 1 в b по х компоненте
	fem.set_boundaries('W', zero );
	fem.set_boundaries('N',
		[&](const Point& p)
		{
			double neighbor_x = b.x - dx;
			if (p.x > neighbor_x)
				return Point{ (p.x - neighbor_x) / dx, 0 };

			return Point{ 0, 0 };
		});
	fem.set_boundaries('E',
		[&](const Point& p)
		{
			double neighbor_y = b.y - dy;
			if (p.y > neighbor_y)
				return Point{ (p.y - neighbor_y) / dy, 0 };

			return Point{ 0, 0 };
		});
	fem.set_boundaries('S', zero);
	
	fem.set_AF(A, F);
	fem.apply_boundaries();
	basis[2 * (n_side_x + n_side_y)] = fem.solve();
	
	// находим базисную функцию для, равную 1 в b по y компоненте
	fem.set_boundaries('W', zero);
	fem.set_boundaries('N',
		[&](const Point& p)
		{
			double neighbor_x = b.x - dx;
		if (p.x > neighbor_x)
			return Point{ 0, (p.x - neighbor_x) / dx};

		return Point{ 0, 0 };
		});
	fem.set_boundaries('E',
		[&](const Point& p)
		{
			double neighbor_y = b.y - dy;
		if (p.y > neighbor_y)
			return Point{ 0, (p.y - neighbor_y) / dy};

		return Point{ 0, 0 };
		});
	fem.set_boundaries('S', zero);

	fem.set_AF(A, F);
	fem.apply_boundaries();
	basis[2 * (n_side_x + n_side_y) + 1] = fem.solve();
	
	//******************ПРАВЫЙ НИЖНИЙ УГОЛ******************
	
	// находим базисную функцию для, равную 1 в {b.x, a.y} по х компоненте
	fem.set_boundaries('W', zero);
	fem.set_boundaries('N', zero);
	fem.set_boundaries('E',
		[&](const Point& p)
		{
			double neighbor_y = a.y + dy;
			if (p.y < neighbor_y)
				return Point{ (neighbor_y - p.y) / dy, 0 };

			return Point{ 0, 0 };
		});
	fem.set_boundaries('S',
		[&](const Point& p)
		{
			double neighbor_x = b.x - dx;
			if (p.x > neighbor_x)
				return Point{ (p.x - neighbor_x) / dx, 0 };

			return Point{ 0, 0 };
		});

	fem.set_AF(A, F);
	fem.apply_boundaries();
	basis[2 * (n_side_x + 2 * n_side_y)] = fem.solve();
	
	// находим базисную функцию для, равную 1 в {b.x, a.y} по y компоненте
	fem.set_boundaries('W', zero);
	fem.set_boundaries('N', zero);
	fem.set_boundaries('E',
		[&](const Point& p)
		{
			double neighbor_y = a.y + dy;
			if (p.y < neighbor_y)
				return Point{ 0, (neighbor_y - p.y) / dy };

			return Point{ 0, 0 };
		});
	fem.set_boundaries('S',
		[&](const Point& p)
		{
			double neighbor_x = b.x - dx;
			if (p.x > neighbor_x)
				return Point{ 0, (p.x - neighbor_x) / dx};

			return Point{ 0, 0 };
		});

	fem.set_AF(A, F);
	fem.apply_boundaries();
	basis[2 * (n_side_x + 2 * n_side_y) + 1] = fem.solve();
	
	//******************ВЕРТИКАЛЬНАЯ СТОРОНА******************
	
	for (int i = 1; i < n_side_y; ++i) {

		//--------------НА ЛЕВОЙ--------------

		// находим базисную функцию для, равную 1 в i по х компоненте
		fem.set_boundaries('W',
			[&](const Point& p)
			{
				double current_y = a.y + i * dy;
				double neighbor_up_y = current_y + dy,
					neighbor_down_y = current_y - dy;
				
				if (neighbor_down_y < p.y && p.y <= current_y )
					return Point{ (p.y - neighbor_down_y) / dy, 0 };

				if ( current_y < p.y && p.y < neighbor_up_y )
					return Point{ (neighbor_up_y - p.y) / dy, 0 };

				return Point{ 0, 0 };
			});
		fem.set_boundaries('N', zero);
		fem.set_boundaries('E', zero);
		fem.set_boundaries('S', zero);

		fem.set_AF(A, F);
		fem.apply_boundaries();
		basis[2 * i] = fem.solve();
		
		// находим базисную функцию для, равную 1 в i по y компоненте
		fem.set_boundaries('W',
			[&](const Point& p)
			{
				double current_y = a.y + i * dy;
				double neighbor_up_y = current_y + dy,
					neighbor_down_y = current_y - dy;

				if (neighbor_down_y < p.y && p.y <= current_y)
					return Point{ 0, (p.y - neighbor_down_y) / dy};

				if (current_y < p.y && p.y < neighbor_up_y)
					return Point{ 0, (neighbor_up_y - p.y) / dy};

				return Point{ 0, 0 };
			});
		fem.set_boundaries('N', zero);
		fem.set_boundaries('E', zero);
		fem.set_boundaries('S', zero);

		fem.set_AF(A, F);
		fem.apply_boundaries();
		basis[2 * i + 1] = fem.solve();
		
		//--------------НА ПРАВОЙ--------------

		// находим базисную функцию для, равную 1 в (n_side_x + n_side_y + i) по х компоненте
		fem.set_boundaries('W', zero);
		fem.set_boundaries('N', zero);
		fem.set_boundaries('E',
			[&](const Point& p)
			{
				double current_y = b.y - i * dy;
				double neighbor_up_y = current_y + dy,
					neighbor_down_y = current_y - dy;

				if (neighbor_down_y < p.y && p.y <= current_y)
					return Point{ (p.y - neighbor_down_y) / dy, 0 };

				if (current_y < p.y && p.y < neighbor_up_y)
					return Point{ (neighbor_up_y - p.y) / dy, 0 };

				return Point{ 0, 0 };
			});
		fem.set_boundaries('S', zero);

		fem.set_AF(A, F);
		fem.apply_boundaries();
		basis[2 * (n_side_x + n_side_y + i)] = fem.solve();
		
		// находим базисную функцию для, равную 1 в i по y компоненте
		fem.set_boundaries('W', zero);
		fem.set_boundaries('N', zero);
		fem.set_boundaries('E',
			[&](const Point& p)
			{
				double current_y = b.y - i * dy;
				double neighbor_up_y = current_y + dy,
					neighbor_down_y = current_y - dy;

				if (neighbor_down_y < p.y && p.y <= current_y)
					return Point{ 0, (p.y - neighbor_down_y) / dy};

				if (current_y < p.y && p.y < neighbor_up_y)
					return Point{ 0, (neighbor_up_y - p.y) / dy };

				return Point{ 0, 0 };
			});
		fem.set_boundaries('S', zero);

		fem.set_AF(A, F);
		fem.apply_boundaries();
		basis[2 * (n_side_x + n_side_y + i) + 1] = fem.solve();
	}

	//******************ГОРИЗОНТАЛЬНАЯ СТОРОНА******************
	
	for (int i = 1; i < n_side_x; ++i) {

		//--------------НА ВЕРХНЕЙ--------------

		// находим базисную функцию для, равную 1 в n_side_x + i по х компоненте
		fem.set_boundaries('W', zero);
		fem.set_boundaries('N',
			[&](const Point& p)
			{
				double current_x = a.x + i * dx;
				double neighbor_right_x = current_x + dx,
					neighbor_left_x = current_x - dx;

				if (neighbor_left_x < p.x && p.x <= current_x)
					return Point{ (p.x - neighbor_left_x) / dx, 0 };

				if (current_x < p.x && p.x < neighbor_right_x)
					return Point{ (neighbor_right_x - p.x) / dx, 0 };

				return Point{ 0, 0 };
			});
		fem.set_boundaries('E', zero);
		fem.set_boundaries('S', zero);

		fem.set_AF(A, F);
		fem.apply_boundaries();
		basis[2 * (n_side_y + i)] = fem.solve();
		
		// находим базисную функцию для, равную 1 в n_side_x + i по y компоненте
		fem.set_boundaries('W', zero);
		fem.set_boundaries('N',
			[&](const Point& p)
			{
				double current_x = a.x + i * dx;
				double neighbor_right_x = current_x + dx,
					neighbor_left_x = current_x - dx;

				if (neighbor_left_x < p.x && p.x <= current_x)
					return Point{ 0, (p.x - neighbor_left_x) / dx};

				if (current_x < p.x && p.x < neighbor_right_x)
					return Point{ 0, (neighbor_right_x - p.x) / dx };

				return Point{ 0, 0 };
			});
		fem.set_boundaries('E', zero);
		fem.set_boundaries('S', zero);

		fem.set_AF(A, F);
		fem.apply_boundaries();
		basis[2 * (n_side_y + i) + 1] = fem.solve();
		
		//--------------НА НИЖНЕЙ--------------

		// находим базисную функцию для, равную 1 в (2 * n_side_y + n_side_x + i) по х компоненте
		fem.set_boundaries('W', zero);
		fem.set_boundaries('N', zero);
		fem.set_boundaries('E', zero);
		fem.set_boundaries('S',
			[&](const Point& p)
			{
				double current_x = b.x - i * dx;
				double neighbor_right_x = current_x + dx,
					neighbor_left_x = current_x - dx;

				if (neighbor_left_x < p.x && p.x <= current_x)
					return Point{ (p.x - neighbor_left_x) / dx, 0 };

				if (current_x < p.x && p.x < neighbor_right_x)
					return Point{ (neighbor_right_x - p.x) / dx, 0 };

				return Point{ 0, 0 };
			});

		fem.set_AF(A, F);
		fem.apply_boundaries();
		basis[2 * (2 * n_side_y + n_side_x + i)] = fem.solve();
		
		// находим базисную функцию для, равную 1 в (2 * n_side_y + n_side_x + i) по y компоненте
		fem.set_boundaries('W', zero);
		fem.set_boundaries('N', zero);
		fem.set_boundaries('E', zero);
		fem.set_boundaries('S',
			[&](const Point& p)
			{
				double current_x = b.x - i * dx;
				double neighbor_right_x = current_x + dx,
					neighbor_left_x = current_x - dx;

				if (neighbor_left_x < p.x && p.x <= current_x)
					return Point{ 0, (p.x - neighbor_left_x) / dx };

				if (current_x < p.x && p.x < neighbor_right_x)
					return Point{ 0, (neighbor_right_x - p.x) / dx };

				return Point{ 0, 0 };
			});

		fem.set_AF(A, F);
		fem.apply_boundaries();
		basis[2 * (2 * n_side_y + n_side_x + i) + 1] = fem.solve();
		
	}

	Matrix W = matrix_form_basis();
	K = W.T().dot(A).dot(W);
}

Matrix FSEM::matrix_form_basis() {
	Matrix W(2 * basis[0].size(), basis.size());
	for (size_t i = 0; i < basis[0].size(); ++i)
		for (size_t j = 0; j < basis.size(); ++j) {
			W[2 * i][j] = basis[j][i].x;
			W[2 * i + 1][j] = basis[j][i].y;
		}
	return W;
}

void FSEM::construct_f_bc2(const std::vector<size_t>& pos,
	const std::vector<vec_function>& g) {
	Matrix W = matrix_form_basis();
	std::vector<double> p_vec(W.size(), 0.0);
	fem.calculate_bc2(pos, g, p_vec);
	f = W.T().dot(p_vec);
}

void FSEM::set_bc1(char side, const vec_function& g) {
	if (side == 'W')
		for (size_t i = 0; i <= n_side_y; ++i) 
			basis_coefficients[i] = g(nodes[i]);

	else if (side == 'N')
		for (size_t i = n_side_y; i <= n_side_x + n_side_y; ++i)
			basis_coefficients[i] = g(nodes[i]);

	else if (side == 'E')
		for (size_t i = n_side_x + n_side_y; i <= n_side_x + 2 * n_side_y; ++i)
			basis_coefficients[i] = g(nodes[i]);

	else 
		for (size_t i = n_side_x + 2 * n_side_y; i < 2 * n_side_x + 2 * n_side_y; ++i)
			basis_coefficients[i] = g(nodes[i]);
}

void FSEM::calculate_coef_Matrix_bc2(const int finish, const int i,
	bool is_cur_Neumann, bool is_prev_Neumann, int& add_N, int& add_D, const int add_basis, Matrix& N,
	Matrix& D) {

	// обработка первой точки на данной границе
	if (is_cur_Neumann && is_prev_Neumann) {
		N[2 * i][add_N] = basis[add_basis][i].x;
		N[2 * i + 1][add_N] = basis[add_basis][i].y;

		N[2 * i][1 + add_N] = basis[1 + add_basis][i].x;
		N[2 * i + 1][1 + add_N] = basis[1 + add_basis][i].y;
	}
	else {
		D[2 * i][add_D] = basis[add_basis][i].x;
		D[2 * i + 1][add_D] = basis[add_basis][i].y;

		D[2 * i][1 + add_D] = basis[1 + add_basis][i].x;
		D[2 * i + 1][1 + add_D] = basis[1 + add_basis][i].y;
	}

	if (is_cur_Neumann && !is_prev_Neumann)
		add_N -= 2;

	for (int j = 2; j < finish; j++) {
		if (is_cur_Neumann) {
			N[2 * i][j + add_N] = basis[j + add_basis][i].x;
			N[2 * i + 1][j + add_N] = basis[j + add_basis][i].y;
		}
		else {
			D[2 * i][j + add_D] = basis[j + add_basis][i].x;
			D[2 * i + 1][j + add_D] = basis[j + add_basis][i].y;
		}
	}

	if (is_cur_Neumann) {
		add_N += finish;
		if (!is_prev_Neumann) {
			add_D += 2;
		}
	}
	else
		add_D += finish;
}

void FSEM::save_bc1(std::vector<double>& coefs_Dirichle, int& dir_id, const int i) {
	coefs_Dirichle[dir_id] = basis_coefficients[i].x;
	coefs_Dirichle[dir_id + 1] = basis_coefficients[i].y;

	dir_id += 2;
}

void FSEM::set_bc2(const std::vector<size_t>& pos,
	const std::vector<vec_function>& g) {

	fem.construct_AF(E, nu, zero);
	auto K_fem  = fem.get_AF().first; // матрица жесткости
	
	size_t n_known_coefs = 0;
	if (!pos[0]) n_known_coefs += n_side_y + 1;
	if (!pos[2]) n_known_coefs += n_side_y + 1;
	if (!pos[1]) n_known_coefs += n_side_x + 1;
	if (!pos[3]) n_known_coefs += n_side_x + 1;

	if (!pos[0] && !pos[1]) n_known_coefs -= 1;
	if (!pos[0] && !pos[3]) n_known_coefs -= 1;
	if (!pos[2] && !pos[1]) n_known_coefs -= 1;
	if (!pos[2] && !pos[3]) n_known_coefs -= 1;

	size_t n_unknown_coefs = 2 * (2 * n_side_x + 2 * n_side_y - n_known_coefs);
	
	// столбцы - значения суперэлементов, соответствующих
	// неизвестным коэффициентам, в узлах мкэ сетки
	Matrix N(K_fem.size(), n_unknown_coefs);

	// столбцы - значения суперэлементов, соответствующих
	// известным коэффициентам, в узлах мкэ сетки
	Matrix D(K_fem.size(), 2 * n_known_coefs);

	// интегралы от ГУ 2 рода * функции формы мкэ
	std::vector<double> p_vec(K_fem.size(), 0.0);

	int finish = 0, prev_pos = 0;
	
	// находим B и C
	for (size_t i = 0; i < fem.psize(); i++) {
		int add_N = 0, add_D = 0, add_basis = 0;

		for (size_t j = 0; j < 4; ++j) {
			if (j == 0)
				prev_pos = 3;
			else
				prev_pos = j - 1;

			if (j == 0 || j == 2)
				finish = 2 * n_side_y;
			else
				finish = 2 * n_side_x;

			calculate_coef_Matrix_bc2(finish, i, pos[j], pos[prev_pos], add_N, add_D, add_basis, N, D);
			add_basis += finish;
		}
	}
	
	fem.calculate_bc2(pos, g, p_vec);
	
	Matrix N_Transposed = N.T();

	Matrix A = N_Transposed.dot(K_fem).dot(N);

	// сохрвняем известные коэффициенты из ГУ Дирихле
	std::vector<double> coefs_Dirichle(2 * n_known_coefs);
	int dir_id = 0;

	if (pos[0] && !pos[3])
		save_bc1(coefs_Dirichle, dir_id, 0);

	if (!pos[0]) {
		for (int i = 0; i < n_side_y; ++i) 
			save_bc1(coefs_Dirichle, dir_id, i);

		if (pos[1]) 
			save_bc1(coefs_Dirichle, dir_id, n_side_y);
	}

	if (!pos[1]) {
		
		for (int i = n_side_y; i < n_side_y + n_side_x; ++i)
			save_bc1(coefs_Dirichle, dir_id, i);

		if (pos[2]) 
			save_bc1(coefs_Dirichle, dir_id, n_side_y + n_side_x);
	}
	
	if (!pos[2]) {
		for (int i = n_side_x + n_side_y; i < n_side_x + 2 * n_side_y; ++i)
			save_bc1(coefs_Dirichle, dir_id, i);

		if (pos[3])
			save_bc1(coefs_Dirichle, dir_id, n_side_x + 2 * n_side_y);
	}

	if (!pos[3]) 
		for (int i = n_side_x + 2 * n_side_y; i < 2 * (n_side_x + n_side_y); ++i) 
			save_bc1(coefs_Dirichle, dir_id, i);
	
	std::vector<double> f_fem;
	if (pos[0] + pos[1] + pos[2] + pos[3] == 4)
		f_fem = N_Transposed.dot(p_vec);
	else
		f_fem = N_Transposed.dot(p_vec) - N_Transposed.dot(K_fem).dot(D).dot(coefs_Dirichle);

	auto [L, U] = LU_decomposition(A);
	std::vector<double> ans = solveLU(L, U, f_fem);
	
	int ans_id = 0;
	if (pos[0])
		for (size_t i = 0; i <= n_side_y; ++i) {
			if (i == 0 && !pos[3])
				continue;
			if (i == n_side_y && !pos[1])
				continue;
			basis_coefficients[i] = { ans[ans_id], ans[ans_id + 1] };
			ans_id += 2;
		}
	
	if (pos[1])
		for (size_t i = n_side_y + 1; i <= n_side_x + n_side_y; ++i) {
			if (i == n_side_x + n_side_y && !pos[2])
				continue;

			basis_coefficients[i] = { ans[ans_id], ans[ans_id + 1] };
			ans_id += 2;
		}

	if (pos[2])
		for (size_t i = n_side_x + n_side_y + 1; i <= n_side_x + 2 * n_side_y; ++i) {
			if (i == n_side_x + 2 * n_side_y && !pos[3])
				continue;
			basis_coefficients[i] = { ans[ans_id], ans[ans_id + 1] };
			ans_id += 2;
		}

	if (pos[3]) {
		for (size_t i = n_side_x + 2 * n_side_y + 1; i < 2 * n_side_x + 2 * n_side_y; ++i) {
			basis_coefficients[i] = { ans[ans_id], ans[ans_id + 1] };
			ans_id += 2;
		}
	}

}

std::vector<Point> FSEM::find_answer() {

	std::vector<Point> res(fem.psize());
	for (size_t i = 0; i < res.size(); i++) 
		for (size_t j = 0; j < nodes.size(); j++) {
			res[i].x += basis_coefficients[j].x * basis[2 * j][i].x +
				basis_coefficients[j].y * basis[2 * j + 1][i].x;
			res[i].y += basis_coefficients[j].x * basis[2 * j][i].y +
				basis_coefficients[j].y * basis[2 * j + 1][i].y;
		}

	return res;
}