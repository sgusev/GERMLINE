// InputManager.cpp: manages input

#include "InputManager.h"
#include "PEDIndividualsExtractor.h"
#include "HMIndividualsExtractor.h"
#include <iostream>
#include <cmath>
using namespace std;


// InputManager(): default constructor
InputManager::InputManager()
{
	// get map of genetic distance
	string choice;
	
	normalPrompt();
	cin >> choice;
	while (!validChoice(choice))
	{
		invalidChoiceMessage();
		normalPrompt();
		cin >> choice;
	}
	setFileFormat(choice);
	instantiatePie();
}

InputManager::~InputManager()
{
	delete pie;
}

bool InputManager::getPhased()
{
	string s; bool phase;
	cout << "Please indicate if the phase is known" << endl
		 << '\t' << "Enter 1 for phased" << endl
		 << '\t' << "Enter 2 for unphased" << endl;
	cin >> s;
	if(s[0] == '1') phase = true;
	else phase = false;

	pie->setPhased(phase);
	return phase;
}

// getIndividuals(): extracts individuals using pie
void InputManager::getIndividuals()
{
	pie->getInput();
}

// normalPrompt():issues prompt asking user which file format
//  will be used for individuals
void InputManager::normalPrompt()
{
	cout << "Please indicate which file format you'll be using for individuals" << endl
		<< "     Enter 1 for PED / Plink" << endl
		<< "     Enter 2 for PHASE / HapMap" << endl;
}

string InputManager::getOutput()
{
	string s;
	cout << "Please indicate output file location: " << endl;
	cin >> s;
	return s;
}

// validChoice(): determines if a string is valid input for file format
bool InputManager::validChoice(string choice)
{
	if (choice[0] == '1' || choice[0] == '2')
		return true;
	else
		return false;
}

PolymorphicIndividualsExtractor* InputManager::getPie()
{
	return pie;
}

// invalidChoiceMessage(): gives message saying format is not valid
void InputManager::invalidChoiceMessage()
{
	cerr << "That choice is not recognized." << endl;
}


// setFileFormat(): sets file format using information from normalPrompt()
void InputManager::setFileFormat(string choice)
{
	if (choice[0] == '1')
		format = PED;
	else if (choice[0] == '2')
		format = HM;
	else
	{
		cerr << "WARNING: InputManager::setFileFormat(): invalid file format choice"
			<< endl;
	}
}


// instantiatePie(): instantiates appropriate pie
void InputManager::instantiatePie()
{
	if (format == PED)
		pie = new PEDIndividualsExtractor();
	else if (format == HM)
		pie = new HMIndividualsExtractor();
}


// end InputManager.cpp
