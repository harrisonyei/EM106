#include "WindowsForm.h"
#include "Matrix.h"
using namespace System;
using namespace System::Windows::Forms;
#ifdef DEBUG_M
int main(){
	Matrix m,m2;
	try{

		Matrix::PrintM(m.Inv());
	} catch(const char* msg){
		std::cout << msg << std::endl;
	}
	system("PAUSE");
	return 0;
}
#else
[STAThread]
void main(array<String^>^ args){
	Application::EnableVisualStyles();
	Application::SetCompatibleTextRenderingDefault(false);
	WindowsFormsApplication_cpp::WindowsForm windowsForm;
	Application::Run(%windowsForm);
}
#endif // DEBUG_M


