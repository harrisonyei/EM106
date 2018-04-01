#include "WindowsForm.h"
#include "Matrix.h"
using namespace System;
using namespace System::Windows::Forms;
#ifdef DEBUG_M
int main(){
	Matrix m;
	try{
		m.Data.push_back({-30,47,47});
		m.Data.push_back({75,-19,6});
		m.Data.push_back({77,7,-18});
		Matrix::Eigen(m,m,m);
	} catch(const char* msg){
		std::cout << msg << std::endl;
	}
	system("PAUSE");
	return 0;
}
#else
[STAThread]
void main(array<String^>^ args)
{
	Application::EnableVisualStyles();
	Application::SetCompatibleTextRenderingDefault(false);
	WindowsFormsApplication_cpp::WindowsForm windowsForm;
	Application::Run(%windowsForm);
}
#endif // DEBUG_M


