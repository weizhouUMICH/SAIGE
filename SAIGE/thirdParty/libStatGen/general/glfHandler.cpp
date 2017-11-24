/*
 *  Copyright (C) 2010  Regents of the University of Michigan
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "glfHandler.h"
#include "BaseQualityHelper.h"

char glfHandler::translateBase[16] = {0, 1, 2, 0, 3, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0};
char glfHandler::backTranslateBase[5] = { 15, 1, 2, 4, 8 };
unsigned char glfHandler::nullLogLikelihoods[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
double glfHandler::nullLikelihoods[10] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};

glfHandler::glfHandler()
{
    isStub = true;
    sections = 0;
    currentSection = 0;
    maxPosition = position = endOfSection = 0;
}

glfHandler::~glfHandler()
{
    // Not safe to close the file here in case a copy of the file was generated
    // if (isOpen())
    //   Close();
}

bool glfHandler::Open(const String & filename)
{
    isStub = false;
    handle = ifopen(filename, "rb");

    if (handle == NULL)
    {
        isStub = true;
        return false;
    }

    if (!ReadHeader())
        ifclose(handle);

    endOfSection = true;

    return handle != NULL;
}

void glfHandler::OpenStub()
{
    isStub = true;
    handle = NULL;

    endOfSection = true;
    data.recordType = 0;
    maxPosition = 1999999999;
    position = maxPosition + 1;
}

bool glfHandler::Create(const String & filename)
{
    isStub = false;
    // glf is in BGZF format.
    handle = ifopen(filename, "wb", InputFile::BGZF);

    if (handle == NULL)
    {
        isStub = true;
        return false;
    }

    WriteHeader();

    return handle != NULL;
}

bool glfHandler::isOpen()
{
    return handle != NULL;
}

bool glfHandler::ReadHeader()
{
    if (isStub)
        return true;

    if (handle == NULL)
        return false;

    char magicNumber[4];

    if (ifread(handle, magicNumber, 4) != 4)
    {
        errorMsg = "unexpected end of file";
        return false;
    }

    if (magicNumber[0] != 'G' || magicNumber[1] != 'L' || magicNumber[2] != 'F')
    {
        errorMsg = "invalid format";
        return false;
    }

    if (magicNumber[3] != 3)
    {
        errorMsg = "unsupported version";
        return false;
    }

    unsigned int headerLength = 0;

    if (ifread(handle, &headerLength, 4) != 4)
    {
        errorMsg = "unexpected end of file";
        return false;
    }

    if (headerLength > 1024 * 1024)
    {
        errorMsg = "header too large -- bailing";
        return false;
    }

    header.SetLength(headerLength + 1);
    header[headerLength] = 0;

    if (headerLength && ifread(handle, header.LockBuffer(headerLength + 1), headerLength) != headerLength)
    {
        errorMsg = "unexpected end of file";
        return false;
    }

    return true;
}

void glfHandler::Close()
{
    if (isOpen())
        ifclose(handle);
}

void glfHandler::Rewind()
{
    if (isOpen())
    {
        ifrewind(handle);

        if (!ReadHeader())
            ifclose(handle);

        endOfSection = true;
    }
}

bool glfHandler::NextSection()
{
    if (isStub)
    {
        endOfSection = true;
        data.recordType = 0;
        maxPosition = 1999999999;
        position = maxPosition + 1;
        return true;
    }

    while (!endOfSection && !ifeof(handle))
        NextEntry();

    endOfSection = false;

    int labelLength = 0;

    currentSection++;
    position = 0;
    if (ifread(handle, &labelLength, sizeof(int)) == sizeof(int))
    {
        ifread(handle, label.LockBuffer(labelLength+1), labelLength * sizeof(char));
        label.UnlockBuffer();

        maxPosition = 0;
        ifread(handle, &maxPosition, sizeof(int));

        return ((maxPosition > 0) && !ifeof(handle));
    }

    return false;
}

bool glfHandler::NextBaseEntry()
{
    bool result = true;

    do
    {
        result = NextEntry();
    }
    while (result && data.recordType == 2);

    return result;
}


bool glfHandler::NextEntry()
{
    if (isStub)
        return false;

    // Read record type
    if (endOfSection || (ifread(handle, &data, 1) != 1))
    {
        endOfSection = true;
        data.recordType = 0;
        position = maxPosition + 1;
        return false;
    }

    // printf("%d/%d\n", data.recordType, data.refBase);

    if (position > maxPosition)
        return true;

    switch (data.recordType)
    {
        case 0 :
            endOfSection = true;
            position = maxPosition + 1;
            return true;
        case 1 :
            if (ifread(handle,((char *) &data) + 1, sizeof(data) - 1) == sizeof(data) - 1)
            {
                data.refBase = translateBase[data.refBase];

                for (int i = 0; i < 10; i++)
                    likelihoods[i] = bQualityConvertor.toDouble(data.lk[i]);

                position = position + data.offset;
                return true;
            }

            // Premature end of file
            data.recordType = 0;
            position = maxPosition + 1;
            return false;
        case 2 :
            while (ifread(handle, ((char *) &data) + 1, sizeof(data) - 4) == sizeof(data) - 4)
            {
                data.refBase = translateBase[data.refBase];

                for (int i = 0; i < 3; i++)
                    likelihoods[i] = bQualityConvertor.toDouble(data.indel.lk[i]);

                position = position + data.offset;

                indelSequence[0].SetLength(abs(data.indel.length[0]) + 1);
                indelSequence[0][abs(data.indel.length[0])] = 0;
                if (ifread(handle, indelSequence[0].LockBuffer(), abs(data.indel.length[0])) != (unsigned int) abs(data.indel.length[0]))
                    break;

                indelSequence[1].SetLength(abs(data.indel.length[1]) + 1);
                indelSequence[1][abs(data.indel.length[1])] = 0;
                if (ifread(handle, indelSequence[1].LockBuffer(), abs(data.indel.length[1])) != (unsigned int) abs(data.indel.length[1]))
                    break;

                return true;
            }

            // Premature end of file
            data.recordType = 0;
            position = maxPosition + 1;
            return false;
    }

    return false;
}

glfEntry & glfEntry::operator = (glfEntry & rhs)
{
    refBase = rhs.refBase;
    recordType = rhs.recordType;
    offset = rhs.offset;
    mapQuality = rhs.mapQuality;

    for (int i = 0; i < 10; i++)
        lk[i] = rhs.lk[i];

    minLLK = rhs.minLLK;
    depth = rhs.depth;

    return * this;
}

const double * glfHandler::GetLikelihoods(int pos)
{
    if (pos == position)
        return likelihoods;

    return nullLikelihoods;
}

const unsigned char * glfHandler::GetLogLikelihoods(int pos)
{
    if (pos == position)
        return data.lk;

    return nullLogLikelihoods;
}

char glfHandler::GetReference(int pos, char defaultBase)
{
    if (pos == position)
        return data.refBase;

    return defaultBase;
}

int glfHandler::GetDepth(int pos)
{
    if (pos == position)
        return data.depth;

    return 0;
}

int glfHandler::GetMapQuality(int pos)
{
    if (pos == position)
        return data.mapQuality;

    return 0;
}

void glfHandler::WriteHeader(const String & headerText)
{
    char magicNumber[4] = {'G', 'L', 'F', 3};

    ifwrite(handle, magicNumber, 4);

    unsigned int headerLength = headerText.Length();

    ifwrite(handle, &headerLength, 4);
    ifwrite(handle, (void *)(const char *) headerText, headerLength);
}

void glfHandler::BeginSection(const String & sectionLabel, int sectionLength)
{
    int labelLength = sectionLabel.Length() + 1;

    ifwrite(handle, &labelLength, sizeof(int));
    ifwrite(handle, (void *)(const char *) sectionLabel, labelLength);
    ifwrite(handle, &sectionLength, sizeof(int));

    label = sectionLabel;
    maxPosition = sectionLength;
}

void glfHandler::EndSection()
{
    char marker = 0;

    ifwrite(handle, &marker, sizeof(char));
}

void glfHandler::WriteEntry(int outputPosition)
{
    data.offset = outputPosition - position;
    position = outputPosition;

    switch (data.recordType)
    {
        case 0 :
            EndSection();
            return;
        case 1 :
            data.refBase = backTranslateBase[data.refBase];
            ifwrite(handle, &data, sizeof(data));
            data.refBase = translateBase[data.refBase];
            return;
        case 2 :
            data.refBase = backTranslateBase[data.refBase];
            ifwrite(handle, &data, sizeof(data) - 3);
            data.refBase = translateBase[data.refBase];

            ifwrite(handle, (void *)(const char *) indelSequence[0], abs(data.indel.length[0]));
            ifwrite(handle, (void *)(const char *) indelSequence[1], abs(data.indel.length[1]));
            return;
    }
}

